#include "simulation.hh"
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <unistd.h>

namespace FCPIC
{
    /********************************************************
    Solving a poisson equation in 2D parallely
    Finite volume discretized equation used is UE + UW + UN + US - 4*UP = f*h^2
     h is the spacing.
     f is the forcing function.
     Nx and Ny are the number of CVs in x and y direction respectively_
    ***********************************************************/

    simulation::simulation(int argc, char **argv)
    {
        aspect = .5; // (INPUT) y_len = aspect (x_len always norm to 1)
        Npart = 500; // (INPUT)
        N = Npart / 10;
        N_int_x = std::sqrt((double)N / aspect);
        N_int_y = aspect * (double)N_int_x;
        N_int_x = 3; //
        N_int_y = 3; //
        N = N_int_x * N_int_y;
        N_x = N_int_x + 2;
        N_y = N_int_y + 2;

        dx = 1 / (double)N_x;
        dy = aspect / (double)N_y;

        MPI_Init(&argc, &argv);

        // retrieve the number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &n_Procs);

        Y_guard_data = new double[N_x];
        X_guard_data = new double[N_int_y];

        Y_guard_data1 = new double[N_x]();
        X_guard_data1 = new double[N_int_y]();
        Y_guard_data2 = new double[N_x]();
        X_guard_data2 = new double[N_int_y]();
        

        bc[X_DIR] = TBD;
        bc[Y_DIR] = TBD;
    }

    simulation::~simulation(){
        delete[] Y_guard_data, X_guard_data,
                 Y_guard_data1, X_guard_data1,
                 Y_guard_data2, X_guard_data2;
    }

    // Creating a virtual cartesian topology
    void simulation::setup_proc_grid()
    {
        // number of processes per row and column
        grid[X_DIR] = 2;
        grid[Y_DIR] = 2;

        if (grid[X_DIR] * grid[Y_DIR] != n_Procs)
            std::cout << "Error MPI: Mismatch of number of processes and process grid" << std::endl;

        int reorder = 1; // reorder process ranks

        // creates a new communicator, grid_comm
        MPI_Cart_create(MPI_COMM_WORLD, 2, grid, wrap_around, reorder, &grid_comm);

        // retrieve new rank and cartesian coordinates of this process
        MPI_Comm_rank(grid_comm, &grid_rank);
        MPI_Cart_coords(grid_comm, grid_rank, 2, grid_coord);

        // calculate ranks of neighboring processes in the grid
        MPI_Cart_shift(grid_comm, 1, 1, &grid_left, &grid_right);
        MPI_Cart_shift(grid_comm, 0, 1, &grid_bottom, &grid_top);

        //Datatype for horizontal data exchange
        MPI_Type_vector(N_int_y, 1, N_x, MPI_DOUBLE, &exchange_field_type[X_DIR]);
        MPI_Type_commit(&exchange_field_type[X_DIR]);

        // Datatype for vertical data exchange
        MPI_Type_vector(N_x, 1, 1, MPI_DOUBLE, &exchange_field_type[Y_DIR]);
        MPI_Type_commit(&exchange_field_type[Y_DIR]);

        //---TESTES---
        //PARA TESTAR COISAS USAR ESTE SÍTIO
        
        std::cout << "Grid rank: " << grid_rank << std::endl;
        std::cout << "N_int_x: " << N_int_x << "   N_int_y: " << N_int_y << std::endl;
        std::cout << "N_x: " << N_x << "   N_y: " << N_y << std::endl;
        double array[25] = { 0, 1, 2, 3, 4,
                             5, 6, 7, 8, 9,
                            10,11,12,13,14,
                            15,16,17,18,19,
                            20,21,22,23,24};
        for(int k = 0; k<25; k++)
            array[k] += grid_rank*25;
        field phi(N_int_x, N_int_y, array);
        //phi.setValue((double) grid_rank);
        sleep(grid_rank);
        phi.print_field(std::cout);
        exchange_charge_buffers(&phi);
        sleep(grid_rank);
        phi.print_field(std::cout);
        
    }

    void simulation::set_periodic_field_bc()
    {
        bc[X_DIR] = PERIODIC;
        bc[Y_DIR] = PERIODIC;
        wrap_around[X_DIR] = 1;
        wrap_around[Y_DIR] = 1;

        setup_proc_grid();
    }

    void simulation::set_conductive_field_bc()
    {
        bc[X_DIR] = CONDUCTIVE;
        bc[Y_DIR] = CONDUCTIVE;
        wrap_around[X_DIR] = 0;
        wrap_around[Y_DIR] = 0;

        setup_proc_grid();
    }

    void simulation::exchange_phi_buffers(field *phi)
    {
        int i, j;

        //Boundary conditions
        if(grid_left == MPI_PROC_NULL){
            if(bc[Y_DIR] == CONDUCTIVE)
            {
                for(j = 0; j<N_int_y; j++)
                    X_guard_data[j] = 0;

                phi->setWestGuard(X_guard_data);
            }
        }

        if(grid_right == MPI_PROC_NULL){
            if(bc[Y_DIR] == CONDUCTIVE)
            {
                for(j = 0; j<N_int_y; j++)
                    X_guard_data[j] = 0;

                phi->setEastGuard(X_guard_data);
            }
        }

        if(grid_top == MPI_PROC_NULL){
            if(bc[X_DIR] == CONDUCTIVE)
            {
                for(int j = 0; j<N_x; j++)
                    Y_guard_data[j] = 0;

                phi->setNorthGuard(Y_guard_data);
            }
        }

        if(grid_bottom == MPI_PROC_NULL){
            if(bc[X_DIR] == CONDUCTIVE)
            {
                for(int j = 0; j<N_x; j++)
                    Y_guard_data[j] = 0;

                phi->setSouthGuard(Y_guard_data);
            }
        }

        i = 1;
        j = 0;

        //Communication stream leftward
        MPI_Sendrecv(&phi->val[WEST_BOUND], 1, exchange_field_type[X_DIR], grid_left, 0, 
                     &phi->val[EAST_GUARD], 1, exchange_field_type[X_DIR], grid_right, 0, 
                     grid_comm, &status);
        
        //Communication stream rightward
        MPI_Sendrecv(&phi->val[EAST_BOUND], 1, exchange_field_type[X_DIR], grid_right, 0, 
                     &phi->val[WEST_GUARD], 1, exchange_field_type[X_DIR], grid_left, 0, 
                     grid_comm, &status);

        //Communication stream upward
        MPI_Sendrecv(&phi->val[NORTH_BOUND], 1, exchange_field_type[Y_DIR], grid_top, 0, 
                     &phi->val[SOUTH_GUARD], 1, exchange_field_type[Y_DIR], grid_bottom, 0, 
                     grid_comm, &status);
        
        //Communication stream downward
        MPI_Sendrecv(&phi->val[SOUTH_BOUND], 1, exchange_field_type[Y_DIR], grid_bottom, 0, 
                     &phi->val[NORTH_GUARD], 1, exchange_field_type[Y_DIR], grid_top, 0, 
                     grid_comm, &status);
    }

    void simulation::exchange_charge_buffers(field *charge)
    {
        int i = 1;
        int j = 0;

        MPI_Sendrecv(&charge->val[NORTH_GUARD], 1, exchange_field_type[Y_DIR], grid_top, 0,  
                     &Y_guard_data[0], N_x, MPI_DOUBLE, grid_bottom, 0, grid_comm, &status);
        if(grid_bottom != MPI_PROC_NULL)
            charge->reduceSouthBound(Y_guard_data);
                     
        MPI_Sendrecv(&charge->val[SOUTH_GUARD], 1, exchange_field_type[Y_DIR], grid_bottom, 0,
                     &Y_guard_data[0], N_x, MPI_DOUBLE, grid_top, 0, grid_comm, &status);
        if(grid_top != MPI_PROC_NULL)
            charge->reduceNorthBound(Y_guard_data);

        MPI_Sendrecv(&charge->val[WEST_GUARD], 1, exchange_field_type[X_DIR], grid_left, 0, 
                     &X_guard_data[0], N_int_y, MPI_DOUBLE, grid_right, 0, grid_comm, &status);
        if(grid_right != MPI_PROC_NULL)
            charge->reduceEastBound(X_guard_data);

        MPI_Sendrecv(&charge->val[EAST_GUARD], 1, exchange_field_type[X_DIR], grid_right, 0, 
                     &X_guard_data[0], N_int_y, MPI_DOUBLE, grid_left, 0, grid_comm, &status);
        if(grid_left != MPI_PROC_NULL)
            charge->reduceWestBound(X_guard_data);
    }

    // Jacobi solver
    void simulation::jacobi(field *phi, field *charge)
    {
        double res, e;
        double global_res = 1.0;
        double tol = 1e-7;

        long int loop = 0;

        // Defining a new temporary field (temp is not part of the domain)
        field temp(N_x, N_y);

        // Starting the iteration loop
        while (global_res > tol)
        {
            // making res 0 so that any error greater than 0 can be equated to this
            res = 0.0;

            // Making the temp field zero after every iteration
            temp.setValue(0.0);

            // exchanges buffer cells
            exchange_phi_buffers(phi);

            for (int i = 1; i <= N_int_y; i++)
                for (int j = 1; j <= N_int_x; j++)
                {
                    temp.val[POSITION] = .25 * (phi->val[NORTH] + phi->val[SOUTH] +
                                                phi->val[EAST] + phi->val[WEST] -
                                                charge->val[POSITION]);

                    e = temp.val[POSITION] - phi->val[POSITION];
                    if (e > res) // norm infty: supremo
                        res = e;
                }

            // Transferring values from temp to u
            for (int i = 1; i <= N_int_y; i++)
                for (int j = 1; j <= N_int_x; j++)
                    phi->val[POSITION] = temp.val[POSITION];

            if (loop % 10 == 0) // balance to be found...
                MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MAX, grid_comm);

            loop++;
        }

        printf("Maximum residual is %e and number of iterations are %ld and I am process %d \n", res, loop, grid_rank);
    }

    void simulation::set_E_value(field *phi, field *Ex_field, field *Ey_field)
    {
        for (int i = 1; i <= N_int_y; i++)
            for (int j = 1; j <= N_int_x; j++)
            {
                Ex_field->val[POSITION] = (phi->val[WEST] - phi->val[EAST]) / (2.f * dx);
                Ey_field->val[POSITION] = (phi->val[SOUTH] - phi->val[WEST]) / (2.f * dy);
            }
        
        if(grid_left == MPI_PROC_NULL){
            if(bc[X_DIR]==CONDUCTIVE){
                for(int i=0; i<=N_y; i++)
                    Ey_field[WEST_GUARD] = 0;
                    Ex_field[WEST_GUARD] = - phi->val[WEST_BOUND] / dx;
            }
        }

        if(grid_right == MPI_PROC_NULL){
            if(bc[X_DIR]==CONDUCTIVE){
                for(int i=0; i<=N_y; i++)
                    Ey_field[EAST_GUARD] = 0;
                    Ex_field[EAST_GUARD] = phi->val[EAST_BOUND] / dx;
            }
        }

        if(grid_top == MPI_PROC_NULL){
            if(bc[Y_DIR]==CONDUCTIVE){
                for(int i=0; i<=N_x; i++)
                    Ex_field[NORTH_GUARD] = 0;
                    Ey_field[NORTH_GUARD] = phi->val[NORTH_BOUND] / dy;
            }
        }

        if(grid_bottom == MPI_PROC_NULL){
            if(bc[Y_DIR]==CONDUCTIVE){
                for(int i=0; i<=N_x; i++)
                    Ex_field[SOUTH_GUARD] = 0;
                    Ey_field[SOUTH_GUARD] = - phi->val[SOUTH_BOUND] / dy;
            }
        }

        exchange_phi_buffers(Ex_field);
        exchange_phi_buffers(Ey_field);
    }

    /*

    // bit dumb: use bool to choose u or charge
    void simulation::write_output_u(domain &subdomain, int rank, int time)
    {
        int l, m;
        int N_local = subdomain.u->N;
        int N_local_x = subdomain.u->Nx;
        int N_local_y = subdomain.u->Ny;
        std::fstream Output_file;

        std::string filename;

        // offset of cell numbering for the subdomain
        offset[X_DIR] = grid_coord[X_DIR] * (N_local_x - 2);
        offset[Y_DIR] = grid_coord[Y_DIR] * (N_local_y - 2);

        filename = "../results/subdomain_" + std::to_string(rank) + "__t_" + std::to_string(time) + ".txt";
        std::string space = "        ";
        Output_file.open(filename, std::ios::out);
        Output_file << "x position" << space << "y position" << space << "field value" << std::endl;

        for (int i = 0; i < N_local; i++)
        {
            if (subdomain.u->bc[i] == NONE)
            {
                l = i % N_local_x;
                m = (int)i / N_local_x;

                double value_x = (l + offset[X_DIR]) * (subdomain.charge->h);
                double value_y = (m + offset[Y_DIR]) * (subdomain.charge->h);

                Output_file << "    " << value_x << space << value_y << space << subdomain.u->val[i] << std::endl;
            }
        }

        Output_file.close();
    }
    void simulation::write_output_charge(domain &subdomain, int rank, int time)
    {
        std::fstream Output_file;
        std::string filename;
        int l, m;
        int N_local = subdomain.u->N;
        int N_local_x = subdomain.u->Nx;
        int N_local_y = subdomain.u->Ny;

        filename = "../results/charge_" + std::to_string(rank) + "__t_" + std::to_string(time) + ".txt";
        std::string space = "   ";

        Output_file.open(filename, std::ios::out);
        int precision = 4;

        for (int i = 0; i < N_local; i++)
        {
            if (i % (N_local_x) == 0)
                Output_file << std::endl;

            Output_file << std::setw(precision) << subdomain.charge->val[i] << space;
            // if (i % 4 == 0)
        }

        Output_file.close();
    }
    */
}