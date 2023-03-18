#include "simulation.hh"
#include "species.hh"

// Inputs: box aspect ratio, box size, no of particles, field BCs, particle BCs

int main(int argc, char **argv)
{
    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::normal_distribution<double> norm(3, 0.5);

    FCPIC::simulation sim(argc, argv);
    sim.set_Xperiodic_field_bc();
    sim.set_Yperiodic_field_bc();

    double t1, t2;
    int l, m;
    t1 = MPI_Wtime();
    /*
    //setup_proc_grid();

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

    // Create Particles Grid
    std::string name = "electron";

    int ppc[2] = {3, 3};
    int range[4] = {0, N_int_local_x, 0, N_int_local_y}; // the latter ones are not included(defined with no guard cells)
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

    std::cout << "****enter loop******" << std::endl;
    //**For Loop Simulation***********************
    for (int counter = 0; counter < 1; counter++)
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

        // Assigning the ghost cells the respective boundary face values for the u field
        set_bc(u);

        // Initial condition for density charge field
        //! print of charge.txt shows err9r: fix
        for (int i = 0; i < N_local; i++)
        {
            if (u->bc[i] == NONE)
            {
                l = i % N_local_x;
                m = (int)i / N_local_x;

                // where all species are added to the charge density grid
                charge->val[i] = test->charge->val[l + m * N_local_x];
            }
            if (u->bc[i] == BUFFER)
            {
                l = i % N_local_x;
                m = (int)i / N_local_x;

                // where all species are added to the charge density grid
                charge->val[i] = test->charge->val[l + m * N_local_x];
            }
        }

        // exchanging buffers of the charge field only once
        // exchange_buffers(charge, N_local_x, N_local_y);

        write_output_charge(subdomain, rank, counter);

        //? Calling jacobi solver
        jacobi(u, N_local_x, N_local_y, charge);

        //? Diagnostics: writes output to the file
        write_output_u(subdomain, rank, counter);
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
    */
    MPI_Finalize();

    return 0;
}