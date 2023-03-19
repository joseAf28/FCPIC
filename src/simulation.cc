#include "simulation.hh"
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>

namespace FCPIC
{
    /********************************************************
    Solving a poisson equation in 2D parallely
    Finite volume discretized equation used is UE + UW + UN + US - 4*UP = f*h^2
     h is the spacing.
     f is the forcing function.
     Nx and Ny are the number of CVs in x and y direction respectively_
    ***********************************************************/


    simulation::simulation(int argc, char **argv){
        aspect = .5; // (INPUT) y_len = aspect (x_len always norm to 1)
        Npart = 500; // (INPUT)
        N= Npart/10;
        N_int_x= std::sqrt((double) N/ aspect);
        N_int_y= aspect*(double) N_int_x;
        N= N_int_x * N_int_y;
        N_x = N_int_x +2;
        N_y = N_int_y +2;

        dx = 1/(double) N_x;
        dy = aspect/(double) N_y;

        MPI_Init(&argc, &argv);

        // retrieve the number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &n_Procs);

        north_recv = new double[N_x];
        north_send = new double[N_x];
        south_recv = new double[N_x];
        south_send = new double[N_x];

        west_recv = new double[N_int_y];
        west_send = new double[N_int_y];
        east_recv = new double[N_int_y];
        east_send = new double[N_int_y];

        bc[X_DIR]= TBD;
        bc[Y_DIR]= TBD;
    }

    simulation::~simulation(){
        delete north_recv, north_send,
               south_recv, south_send,
               west_recv, west_send,
               east_recv, east_send;
    }
    
    // Creating a virtual cartesian topology
    void simulation::setup_proc_grid()
    {
        // number of processes per row and column
        grid[X_DIR] = 2;
        grid[Y_DIR] = 2;

        if (grid[X_DIR] * grid[Y_DIR] != n_Procs)
            std::cout << "Error MPI: Mismatch of number of processes and process grid" << std::endl;

        int reorder = 1;            // reorder process ranks

        // creates a new communicator, grid_comm
        MPI_Cart_create(MPI_COMM_WORLD, 2, grid, wrap_around, reorder, &grid_comm);

        // retrieve new rank and cartesian coordinates of this process
        MPI_Comm_rank(grid_comm, &grid_rank);
        MPI_Cart_coords(grid_comm, grid_rank, 2, grid_coord);

        // calculate ranks of neighboring processes in the grid
        MPI_Cart_shift(grid_comm, 1, 1, &grid_left, &grid_right);
        MPI_Cart_shift(grid_comm, 0, 1, &grid_bottom, &grid_top);

        //---TESTES---
        //PARA TESTAR COISAS USAR ESTE SÍTIO
        /*
        std::cout << "Grid rank: " << grid_rank << std::endl;
        std::cout << "N_int_x: " << N_int_x << "   N_int_y: " << N_int_y << std::endl;
        std::cout << "N_x: " << N_x << "   N_y: " << N_y << std::endl;
        field phi(N_int_x, N_int_y);
        phi.setValue((double) grid_rank);
        exchange_phi_buffers(&phi);
        phi.print_field(std::cout);
        if(grid_rank == -1){
            phi.print_field(std::cout);
        }
        */

    }

    void simulation::set_Xperiodic_field_bc()
    {
        bc[X_DIR] = PERIODIC;
        wrap_around[X_DIR] = 1;

        if(bc[Y_DIR] != TBD)
            setup_proc_grid();
    }

    void simulation::set_Yperiodic_field_bc()
    {
        bc[Y_DIR] = PERIODIC;
        wrap_around[Y_DIR] = 1;

        if(bc[X_DIR] != TBD)
            setup_proc_grid();
    }

    void simulation::set_Xconductive_field_bc()
    {
        bc[X_DIR] = CONDUCTIVE;
        wrap_around[X_DIR] = 0;

        if(bc[Y_DIR] != TBD)
            setup_proc_grid();
    }

    void simulation::set_Yconductive_field_bc()
    {
        bc[Y_DIR] = CONDUCTIVE;
        wrap_around[Y_DIR] = 0;

        if(bc[X_DIR] != TBD)
            setup_proc_grid();
    }

    // communication between the processes
    void simulation::exchange_phi_buffers(field *phi)
    {
        if(grid_left != MPI_PROC_NULL){
            phi->getWestBound(west_send);

            MPI_Sendrecv(west_send, N_int_y, MPI_DOUBLE, grid_left, 0, 
                        west_recv, N_int_y, MPI_DOUBLE, grid_left, 0, 
                        grid_comm, &status);

            phi->setWestGuard(west_recv);
        }
        else{
            if(bc[Y_DIR] == CONDUCTIVE)
            {
                for(int j = 0; j<N_int_y; j++)
                    west_send[j] = 0;

                phi->setWestGuard(west_send);
            }
        }

        if(grid_right != MPI_PROC_NULL){
            phi->getEastBound(east_send);

            MPI_Sendrecv(east_send, N_int_y, MPI_DOUBLE, grid_right, 0, 
                        east_recv, N_int_y, MPI_DOUBLE, grid_right, 0, 
                        grid_comm, &status);

            phi->setEastGuard(east_recv);
        }
        else{
            if(bc[Y_DIR] == CONDUCTIVE)
            {
                for(int j = 0; j<N_int_y; j++)
                    east_send[j] = 0;

                phi->setEastGuard(east_send);
            }
        }

        if(grid_top != MPI_PROC_NULL){
            phi->getNorthBound(north_send);

            MPI_Sendrecv(north_send, N_x, MPI_DOUBLE, grid_top, 0, 
                        north_recv, N_x, MPI_DOUBLE, grid_top, 0, 
                        grid_comm, &status);

            phi->setNorthGuard(north_recv);
        }
        else{
            if(bc[X_DIR] == CONDUCTIVE)
            {
                for(int j = 0; j<N_x; j++)
                    north_send[j] = 0;

                phi->setNorthGuard(north_send);
            }
        }

        if(grid_bottom != MPI_PROC_NULL){
            phi->getSouthBound(south_send);

            MPI_Sendrecv(south_send, N_x, MPI_DOUBLE, grid_bottom, 0, 
                        south_recv, N_x, MPI_DOUBLE, grid_bottom, 0, 
                        grid_comm, &status);

            phi->setSouthGuard(south_recv);
        }
        else{
            if(bc[X_DIR] == CONDUCTIVE)
            {
                for(int j = 0; j<N_x; j++)
                    south_send[j] = 0;

                phi->setSouthGuard(south_send);
            }
        }
    }

    // Jacobi solver
    void simulation::jacobi(field *phi, field *charge)
    {
        double res, e;
        int l, m;
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

            for(int i=1; i<=N_int_y; i++)
                for(int j=1; j<=N_int_x; j++)
                {
                    temp.val[POSITION] = .25*(phi->val[NORTH] + phi->val[SOUTH] +
                                              phi->val[EAST] + phi->val[WEST] -
                                              charge->val[POSITION]);
                    
                    e = temp.val[POSITION] - phi->val[POSITION];
                    if (e > res) // norm infty: supremo
                        res = e;
                }
            
            // Transferring values from temp to u
            for(int i=1; i<=N_int_y; i++)
                for(int j=1; j<=N_int_x; j++)
                    phi->val[POSITION] = temp.val[POSITION];


            if (loop % 10 == 0) // balance to be found...
                MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MAX, grid_comm);

            loop++;
        }

        printf("Maximum residual is %e and number of iterations are %ld and I am process %d \n", res, loop, grid_rank);
    }

    void simulation::set_E_value(field *phi, field *Ex_field, field *Ey_field){
        for(int i=1; i<=N_int_y; i++)
            for(int j=1; j<=N_int_x; j++){
                Ex_field->val[POSITION] = (phi->val[WEST]-phi->val[EAST])/(2.f*dx);
                Ey_field->val[POSITION] = (phi->val[NORTH]-phi->val[SOUTH])/(2.f*dy);
            }
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