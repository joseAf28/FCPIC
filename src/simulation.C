#include "simulation.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

namespace simulation
{
    /********************************************************
    Solving a poisson equation in 2D parallely
    Finite volume discretized equation used is UE + UW + UN + US - 4*UP = f*h^2
     h is the spacing.
     f is the forcing function.
     Nx and Ny are the number of CVs in x and y direction respectively_
    ***********************************************************/

    // MPI variables
    int rank;                                         // rank of the current process
    int grid_rank;                                    // rank of the current proces in the virtual grid
    int grid_top, grid_bottom, grid_left, grid_right; // ranks of the neighbouring processes
    int n_Procs;                                      // total number of processes
    int grid[2];                                      // virtual grid dimensions
    MPI_Datatype exchange_buffer_type[2];             // MPI_datatype for exchange of buffer cell data
    MPI_Comm grid_comm;                               // grid COMMUNICATOR
    int offset[2];                                    // offset for cell numbering for subdomains
    int grid_coord[2];                                // coordinates of the process in the virtual grid
    MPI_Status status;

    // Creating a virtual cartesian topology
    void setup_proc_grid()
    {
        int wrap_around[2];
        int reorder;

        // retrieve the number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &n_Procs);

        // number of processes per row and column
        grid[X_DIR] = 2;
        grid[Y_DIR] = 2;

        if (grid[X_DIR] * grid[Y_DIR] != n_Procs)
            std::cout << "Error MPI: Mismatch of number of processes and process grid" << std::endl;

        // creating a virtual process grid topology
        wrap_around[X_DIR] = 0;
        wrap_around[Y_DIR] = 0; // dont connect first and last processes
        reorder = 1;            // reorder process ranks

        // creates a new communicator, grid_comm
        MPI_Cart_create(MPI_COMM_WORLD, 2, grid, wrap_around, reorder, &grid_comm);

        // retrieve new rank and cartesian coordinates of this process
        MPI_Comm_rank(grid_comm, &grid_rank);
        MPI_Cart_coords(grid_comm, grid_rank, 2, grid_coord);

        // calculate ranks of neighboring processes in the grid
        MPI_Cart_shift(grid_comm, 1, 1, &grid_left, &grid_right);
        MPI_Cart_shift(grid_comm, 0, 1, &grid_bottom, &grid_top);
    }

    // Packing data for communication across processes
    void setup_MPI_datatypes(int N_int_x, int N_int_y, int N_x)
    {
        // Datatype for horizontal data exchange
        MPI_Type_vector(N_int_y, 1, N_x, MPI_DOUBLE, &exchange_buffer_type[X_DIR]);
        MPI_Type_commit(&exchange_buffer_type[X_DIR]);

        // Datatype for vertical data exchange
        MPI_Type_vector(N_int_x, 1, 1, MPI_DOUBLE, &exchange_buffer_type[Y_DIR]);
        MPI_Type_commit(&exchange_buffer_type[Y_DIR]);
    }

    //! very important - Setting flags for ghost and buffer cells
    void set_ghost_buffer_flag(domain &subdomain)
    {
        int i, l, m;

        // to simplify notation
        field *u = subdomain.u;

        int N_Cells = u->Nx * u->Ny;
        int N_Cells_x = u->Nx;
        int N_Cells_y = u->Ny;

        for (i = 0; i < N_Cells; i++) // identify in which BC we are: dirichelet, buffer, none(bulk domain)
        {
            l = i % N_Cells_x;      // gives the position of the CV along x
            m = (int)i / N_Cells_x; // gives the position of the CV along y

            if ((l == 0 && grid_left == MPI_PROC_NULL) || (l == N_Cells_x - 1 && grid_right == MPI_PROC_NULL) || (m == 0 && grid_bottom == MPI_PROC_NULL) || (m == N_Cells_y - 1 && grid_top == MPI_PROC_NULL))
            {
                u->bc[i] = PHYSICAL;
            }
            else if (l == 0 || l == N_Cells_x - 1 || m == 0 || m == N_Cells_y - 1)
            {
                u->bc[i] = BUFFER;
            }
            else
            {
                u->bc[i] = NONE;
            }
        }
    }

    // Setting boundary condition values to boundary cells
    void set_bc(field *phi)
    {

        int l, m;
        int N_Cells_x = phi->Nx;
        int N_Cells_y = phi->Ny;
        int N_Cells = N_Cells_x * N_Cells_y;

        for (int i = 0; i < N_Cells; i++)
        {
            l = i % N_Cells_x;
            m = (int)i / N_Cells_x;

            if (phi->bc[i] == PHYSICAL)
            {
                if (l == 0)
                    phi->val[i] = phi->bc_val[XMIN];
                else if (l == N_Cells_x - 1)
                    phi->val[i] = phi->bc_val[XMAX];
                else if (m == 0)
                    phi->val[i] = phi->bc_val[YMIN];
                else if (m == N_Cells_y - 1)
                    phi->val[i] = phi->bc_val[YMAX];
            }
        }
    }

    // communication between the processes
    void exchange_buffers(field *phi, int N_local_x, int N_local_y)
    {
        // All traffic in direction "top"
        MPI_Sendrecv(&(phi->val[(N_local_y - 2) * (N_local_x) + 1]), 1, exchange_buffer_type[Y_DIR], grid_top, 0, &(phi->val[1]), 1, exchange_buffer_type[Y_DIR], grid_bottom, 0, grid_comm, &status);

        // All traffic in direction "bottom"
        MPI_Sendrecv(&(phi->val[N_local_x + 1]), 1, exchange_buffer_type[Y_DIR], grid_bottom, 0, &(phi->val[(N_local_y - 1) * (N_local_x) + 1]), 1, exchange_buffer_type[Y_DIR], grid_top, 0, grid_comm, &status);

        // All traffic in direction "left"
        MPI_Sendrecv(&(phi->val[N_local_x + 1]), 1, exchange_buffer_type[X_DIR], grid_left, 0, &(phi->val[(2 * N_local_x) - 1]), 1, exchange_buffer_type[X_DIR], grid_right, 0, grid_comm, &status);

        // All traffic in direction "right"
        MPI_Sendrecv(&(phi->val[2 * (N_local_x - 1)]), 1, exchange_buffer_type[X_DIR], grid_right, 0, &(phi->val[N_local_x]), 1, exchange_buffer_type[X_DIR], grid_left, 0, grid_comm, &status);
    }

    // Jacobi solver
    void jacobi(field *phi, int Nx, int Ny, field *charge)
    {
        double res, e;
        int l, m;
        int N = Nx * Ny;
        double global_res = 1.0;
        double tol = 1e-7;

        long int loop = 0;

        // Defining a new temporary field (temp is not part of the domain)
        field temp(Nx, Ny);

        // Starting the iteration loop
        while (global_res > tol)
        {
            // making res 0 so that any error greater than 0 can be equated to this
            res = 0.0;

            // Making the temp field zero after every iteration
            temp.set_field_value(0.0);

            // exchanges buffer cells
            exchange_buffers(phi, Nx, Ny);

            double u_E, u_W, u_N, u_S, u_P;

            for (int i = 0; i < N; i++)
            {
                if (phi->bc[i] == NONE)
                {
                    l = i % Nx;
                    m = (int)i / Nx;

                    u_E = phi->val[RIGHT];
                    u_W = phi->val[LEFT];
                    u_N = phi->val[UP];
                    u_S = phi->val[DOWN];
                    u_P = phi->val[P];

                    if (l == 1 && grid_left == MPI_PROC_NULL)
                        temp.val[P] += 2.0 * u_W - u_P;
                    else
                        temp.val[P] += u_W;
                    if (l == Nx - 2 && grid_right == MPI_PROC_NULL)
                        temp.val[P] += 2.0 * u_E - u_P;
                    else
                        temp.val[P] += u_E;
                    if (m == 1 && grid_bottom == MPI_PROC_NULL)
                        temp.val[P] += 2.0 * u_S - u_P;
                    else
                        temp.val[P] += u_S;
                    if (m == Ny - 2 && grid_top == MPI_PROC_NULL)
                        temp.val[P] += 2.0 * u_N - u_P;
                    else
                        temp.val[P] += u_N;

                    temp.val[P] -= charge->val[P] * charge->h * charge->h;
                    temp.val[P] = temp.val[P] / 4.0;

                    e = temp.val[P] - phi->val[P];
                    if (e > res) // norm infty: supremo
                        res = e;
                }
            }

            // Transferring values from temp to u
            for (int i = 0; i < N; i++)
            {
                if (phi->bc[i] == NONE)
                    phi->val[i] = temp.val[i];
            }

            if (loop % 10 == 0) // balance to be found...
                MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MAX, grid_comm);

            loop++;
        }

        printf("Maximum residual is %e and number of iterations are %ld and I am process %d \n", res, loop, rank);
    }

    void write_output(domain &subdomain, int rank, int time)
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

        filename = "../outputs/subdomain_" + std::to_string(rank) + "__t_" + std::to_string(time) + ".txt";
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
}