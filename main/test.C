#include "simulation.h"
#include "species.h"

using namespace simulation;

int main(int argc, char **argv)
{
    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::normal_distribution<double> norm(3, 0.5);

    MPI_Init(&argc, &argv);

    double t1, t2;
    int l, m;
    t1 = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    setup_proc_grid();

    // Box Simulation Info
    //***********************************
    //  Defining the number of total interior control volumes
    int N_int_x = 100;
    int N_int_y = 100;

    // With ghost cells
    int N_Cells_x = N_int_x + 2;
    int N_Cells_y = N_int_y + 2;
    int N_Cells = N_Cells_x * N_Cells_y;

    // Defining the number of interior local control volumes
    int N_int_local_x = N_int_x / grid[X_DIR];
    int N_int_local_y = N_int_y / grid[Y_DIR];

    // With buffer cells
    int N_local_x = N_int_local_x + 2;
    int N_local_y = N_int_local_y + 2;
    int N_local = N_local_x * N_local_y;
    //*************************************

    // Setting up MPI datatypes for exchange of values in buffer in x and y direction.
    setup_MPI_datatypes(N_int_local_x, N_int_local_y, N_local_x);

    // Create Particles Grid
    std::string name = "electron";

    int ppc[2] = {3, 3};
    int range[4] = {0, N_int_local_x, 0, N_int_local_y}; // the latter ones are not included
    int vec_u[2] = {0, 20};

    float vf[3] = {0.1, 0.1, 0.1};
    float vth[3] = {0.1, 0.1, 0.1};

    species *test = new species(name, ppc, range, vec_u, vf, vth);
    test->set_X();
    test->set_U();
    test->get_charge();

    // test.write_output(proc_rank);

    // Defining a domain of type Domain
    field *u = new field(N_local_x, N_local_y);
    field *charge = new field(N_local_x, N_local_y);

    // subdomain: to organize the objects: links th Potential Field and the Charge Density Field
    domain subdomain(u, charge);

    // Assigning different flags to different boundary conditions
    set_ghost_buffer_flag(subdomain);

    std::cout << "****enter loop******" << std::endl;
    //**For Loop Simulation***********************
    for (int counter = 0; counter < 2; counter++)
    {
        std::cout << "counter: " << counter << std::endl;

        //? get density charge field update
        // test.get_charge();

        // Initializing all boundary condition values to zero: in or out of the loop: idk yet
        for (int i = 0; i < 5; i++)
            u->bc_val[i] = 0.0;

        // Setting up boundary condition values: change later

        u->bc_val[XMAX] = 10.0;
        u->bc_val[XMIN] = 10.0;
        u->bc_val[YMAX] = 10.0;
        u->bc_val[YMIN] = 10.0;

        // initialization of local CVs in the process
        u->set_field_value(0.0);

        // Assigning the ghost cells the respective boundary face values
        set_bc(u);

        // Initial condition for density charge field
        for (int i = 0; i < N_local; i++)
        {
            if (u->bc[i] == NONE || grid_left == MPI_PROC_NULL || grid_right == MPI_PROC_NULL || grid_top == MPI_PROC_NULL || grid_bottom == MPI_PROC_NULL)
            {
                l = i % N_local_x;
                m = (int)i / N_local_x;

                // where all species are added to the chrage density grid
                if (counter == 0)
                    charge->val[i] = test->charge->val[i];

                // simulate the change in the density profile as a result of the particle pusher
                if (counter == 1)
                {
                    charge->val[i] = test->charge->val[i];
                }
            }
        }

        // exchanging buffers of the charge field only once
        exchange_buffers(charge, N_local_x, N_local_y);

        //? Calling jacobi solver
        jacobi(u, N_local_x, N_local_y, charge);

        //? Diagnostics: writes output to the file
        write_output(subdomain, rank, counter);
        test->write_output(rank, counter);

        //? Field Interpolation
        //? Particle Pusher - Change in Particle Grid
        for (int i = 0; i < test->charge->val.size(); i++)
        {
            test->charge->val[i] = norm(rng);
        }
    }
    //*********end loop simulation**********

    t2 = MPI_Wtime();

    printf("Total time taken is %g\n", t2 - t1);

    delete u;
    delete charge;
    delete test;
    MPI_Finalize();

    return 0;
}