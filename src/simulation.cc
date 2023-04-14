// FCPIC - 2D Particle-in-Cell code using MPI
// Guilherme Crispim, João Palma, José Afonso
// Advanced Topics in Computational Physics, 2023, IST

// File simulation.cc:
// Implementation of constructor, destructor and simulation numerical
// methods inside the class Simulation

#include "simulation.hh"
#include <fstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <unistd.h>

namespace FCPIC
{
    /*
        Constructor of class Simulation
        Inputs: int main argument number and list
        + Start MPI environment
        + Call respective functions for getting simulation parameters,
        setting the spatial grid dimensions, and setting up the MPI
        communication grid
    */
    simulation::simulation(int argc, char **argv) : FCPIC_base(),
                                                    total_time(0.), setup_time(0.), hdf5_time(0.), particle_time(0.), field_time(0.)
    {
        MPI_Init(&argc, &argv);
        setTime(); // Time count

        MPI_Comm_size(MPI_COMM_WORLD, &n_Procs); // Number of processes
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);    // First rank attribution
        readArgs(argc, argv);                    // Read user inputs
        setParams();                             // Set all relevant parameters

        if (rank == 0)
            printTitle();

        // Aux arrays for MPI grid communication (charge and electric potential)
        Y_guard_data = new float[N_x];
        X_guard_data = new float[N_int_y];

        setup_proc_grid();   // Set the MPI communication into 2D process grid
        setTime(setup_time); // Count the elapsed time to "setup"
        confirmParams();     // User confirmation of params
        setTime();           // Discard time elapsed in user interaction
    }

    /*
        Destructor of class Simulation
        + Free all memory used in intermediate arrays
    */
    simulation::~simulation()
    {
        // Free all memory allocated
        MPI_Type_free(&exchange_field_type[X_DIR]);
        MPI_Type_free(&exchange_field_type[Y_DIR]);
        MPI_Type_free(&exchange_part_type);

        MPI_Finalize();

        delete[] Y_guard_data, X_guard_data;
    }

    /*
        Function get_charge
        Inputs: charge field object and species object
        + Deposits charge from a single species into the charge
        discrete grid. This is done using linear interpolation
    */
    void simulation::get_charge(field *charge, species *part)
    {
        // Aux variables
        int i, j;
        float wx, wy;
        // Equilibrium density, for normalization
        float density0 = ((float)part->np_sim) / ((float)N_total_x * (float)N_total_y);
        float q = part->q; // Charge
        for (int k = 0; k < part->np; k++)
        {
            i = part->vec[k].iy;
            j = part->vec[k].ix;
            wx = part->vec[k].x;
            wy = part->vec[k].y;

            // Depositing charge by linear interpolation, from the areas of rectangles
            // described by the particle position and the 4 adjacent cell corners
            charge->val[POSITION] += (dx - wx) * (dy - wy) * q / (dx * dy * density0);
            charge->val[EAST] += wx * (dy - wy) * q / (dx * dy * density0);
            charge->val[NORTH] += (dx - wx) * wy * q / (dx * dy * density0);
            charge->val[NORTHEAST] += wx * wy * q / (dx * dy * density0);
        }

        // Including the effect of a static "ionic" background
        for (int i = 1; i <= N_int_y; i++)
            for (int j = 1; j <= N_int_x; j++)
                charge->val[POSITION] -= q;
    }

    /*
        Function jacobi
        Inputs: charge and electric potential field objects
        + Uses the Jacobi iterative method to solve the Poisson equation
        and get the electric potential from the charge density
    */
    void simulation::jacobi(field *phi, field *charge)
    {
        float res, e;           // Local error (single process)
        float global_res = 1.0; // Global error (between all processes)
        float tol = 1e-4;       // Tolerance of the method

        long int loop = 0;

        phi->setValue(0.0);
        // Defining a new temporary field (temp is not part of the domain)
        field temp(this);

        // Starting the iteration loop
        while (global_res > tol)
        {

            // Making res 0 so that any error greater than 0 can be equated to this
            res = 0.0;

            // Making the temp field zero after every iteration
            temp.setValue(0.0);

            // Exchanges buffer cells to include information of the adjacent processes
            exchange_phi_buffers(phi);

            float dx2 = dx*dx;
            float dy2 = dy*dy;
            float d2 = dx2*dy2;
            float invd2 = .5/(dx2+dy2);

            for (int i = 1; i <= N_int_y; i++)
                for (int j = 1; j <= N_int_x; j++)
                {
                    // In the case of periodic conditions, the system is not totally defined:
                    // the BCs do not fix an absolute value, only fix relative values between
                    // boundaries. Therefore, we choose the bottom left point to be 0
                    if (i == 1 && j == 1 && grid_rank == 0 && bc == PERIODIC)
                        temp.val[POSITION] = 0;
                    // Application of the Jacobi iteration formula
                    else
                        //temp.val[POSITION] = .25 * (phi->val[NORTH] + phi->val[SOUTH] + phi->val[EAST] + phi->val[WEST] + dx * dx * charge->val[POSITION]);
                        temp.val[POSITION] = invd2 * (dy2*(phi->val[NORTH] + phi->val[SOUTH]) + dx2*(phi->val[EAST] + phi->val[WEST]) + d2 * charge->val[POSITION]);

                    // Getting local max deviation between iterations
                    e = fabs(temp.val[POSITION] - phi->val[POSITION]);
                    if (e > res)
                        res = e;
                }

            // Transferring values from temp to u
            for (int i = 1; i <= N_int_y; i++)
                for (int j = 1; j <= N_int_x; j++)
                    phi->val[POSITION] = temp.val[POSITION];

            // Communication to obtain the global max deviation
            if (loop % 10 == 0) // balance to be found...
                MPI_Allreduce(&res, &global_res, 1, MPI_FLOAT, MPI_MAX, grid_comm);

            loop++;
        }

        exchange_phi_buffers(phi);
    }

    /*
        Function set_E_value
        Inputs: electric potential and electric field field objects
        + Calculate the electric field grid from the electric potential
    */
    void simulation::set_E_value(field *phi, field *Ex_field, field *Ey_field)
    {
        // Central finite difference for internal points
        for (int i = 1; i <= N_int_y; i++)
            for (int j = 1; j <= N_int_x; j++)
            {
                Ex_field->val[POSITION] = (phi->val[WEST] - phi->val[EAST]) / (2.f * dx);
                Ey_field->val[POSITION] = (phi->val[SOUTH] - phi->val[NORTH]) / (2.f * dy);
            }

        // For conductive boundary conditions, the E-field must be perpendicular to the surface.
        // Parallel components are set to 0, perpendicular are defined from forward/backward
        // finite differences
        if (grid_left == MPI_PROC_NULL)
        {
            if (bc == CONDUCTIVE)
            {
                for (int i = 0; i < N_y; i++)
                {
                    Ey_field->val[WEST_GUARD] = 0;
                    Ex_field->val[WEST_GUARD] = -phi->val[WEST_BOUND] / dx;
                }
            }
        }

        if (grid_right == MPI_PROC_NULL)
        {
            if (bc == CONDUCTIVE)
            {
                for (int i = 0; i < N_y; i++)
                {
                    Ey_field->val[EAST_GUARD] = 0;
                    Ex_field->val[EAST_GUARD] = phi->val[EAST_BOUND] / dx;
                }
            }
        }

        if (grid_top == MPI_PROC_NULL)
        {
            if (bc == CONDUCTIVE)
            {
                for (int j = 0; j < N_x; j++)
                {
                    Ex_field->val[NORTH_GUARD] = 0;
                    Ey_field->val[NORTH_GUARD] = phi->val[NORTH_BOUND] / dy;
                }
            }
        }

        if (grid_bottom == MPI_PROC_NULL)
        {
            if (bc == CONDUCTIVE)
            {
                for (int j = 0; j < N_x; j++)
                {
                    Ex_field->val[SOUTH_GUARD] = 0;
                    Ey_field->val[SOUTH_GUARD] = -phi->val[SOUTH_BOUND] / dy;
                }
            }
        }

        exchange_phi_buffers(Ex_field);
        exchange_phi_buffers(Ey_field);
    }

    /*
        Function field_interpolate
        Inputs: electric field objects, particle struct and two float references for
        writing the interpolated fields
        + Interpolate the electric field grid to obtain the values at the particle position
    */
    void simulation::field_interpolate(field *Ex, field *Ey, float &Ex_i, float &Ey_i, part *prt)
    {
        int i = prt->iy;
        int j = prt->ix;

        float wx = prt->x;
        float wy = prt->y;

        // Areas of the rectangles ("undoing" the process of get_charge)
        float A_pos = (dx - wx) * (dy - wy);
        float A_e = wx * (dy - wy);
        float A_n = (dx - wx) * wy;
        float A_ne = wx * wy;

        // Weighting the values with the areas
        Ex_i = A_pos * Ex->val[POSITION] +
               A_e * Ex->val[EAST] +
               A_n * Ex->val[NORTH] +
               A_ne * Ex->val[NORTHEAST];

        Ey_i = A_pos * Ey->val[POSITION] +
               A_e * Ey->val[EAST] +
               A_n * Ey->val[NORTH] +
               A_ne * Ey->val[NORTHEAST];

        Ex_i /= (dx * dy);
        Ey_i /= (dx * dy);
    }

    /*
        Function init_pusher
        Inputs: electric field objects, species object
        + Push the velocity half a step backwards to stagger them in reference
        to the positions (as needed for the BORIS pusher)
    */
    void simulation::init_pusher(field *Ex, field *Ey, species *spec)
    {

        int Np = spec->np;
        float q = spec->q;
        float m = spec->m;
        for (int i = 0; i < Np; i++)
        {
            float Ex_i = 0.f;
            float Ey_i = 0.f;

            // Get the electric field at each particle position
            field_interpolate(Ex, Ey, Ex_i, Ey_i, &spec->vec[i]);

            // Use the velocity advancing equation, replacing dt->-dt/2
            spec->vec[i].ux += -0.5 * q / m * Ex_i * dt;
            spec->vec[i].uy += -0.5 * q / m * Ey_i * dt;
        }
    }

    /*
        Function particle_pusher
        Inputs: electric field objects, species object
        + Push the particle velocities and positions according to the
        BORIS method
    */
    void simulation::particle_pusher(field *Ex, field *Ey, species *spec)
    {

        int Np = spec->np;
        float q = spec->q;
        float m = spec->m;
        for (int i = 0; i < Np; i++)
        {
            float Ex_i = 0.f;
            float Ey_i = 0.f;

            // Get the electric field at each particle position
            field_interpolate(Ex, Ey, Ex_i, Ey_i, &spec->vec[i]);

            // Advance the velocities
            spec->vec[i].ux += q / m * Ex_i * dt;
            spec->vec[i].uy += q / m * Ey_i * dt;

            // Advance the positions
            spec->vec[i].x += spec->vec[i].ux * dt;
            spec->vec[i].y += spec->vec[i].uy * dt;
        }
    }

    /*
        Function setTime
        Inputs: time variable reference
        + Add, to the reference given, the elapsed
        time since last checking
    */
    void simulation::setTime(float &task_time)
    {
        time2 = MPI_Wtime();
        task_time += time2 - time1;
        total_time += time2 - time1;
        time1 = time2;
    }

    /*
        Function setTime
        + Ignore away elapsed time since last checking
    */
    void simulation::setTime()
    {
        time1 = MPI_Wtime();
    }

    /*
        Function runSimulation
        Inputs: electric field, potential and charge objects,
        reference to vector of species and file name for writing
        outputs
        + Runs sequencially all methods of the simulation throughout
        its entire time, evoking as well all output writing functions
        (HDF5), time monitoring and IO.
    */
    void simulation::run_simulation(field *Ex, field *Ey,
                                    field *phi, field *charge,
                                    std::vector<species> &spec_vec,
                                    std::string filename)
    {
        setupHDF5(filename); // Prepare HDF5 writing

        setTime(setup_time);

        // Charge acquisition
        for (int i = 0; i < Nspecies; i++)
            get_charge(charge, &spec_vec[i]);

        exchange_charge_buffers(charge);
        jacobi(phi, charge);
        set_E_value(phi, Ex, Ey);

        setTime(field_time);

        for (int i = 0; i < Nspecies; i++)
            init_pusher(Ex, Ey, &spec_vec[i]);

        setTime(particle_time);

        for (int counter = 0; 1; counter++)
        {
            // Testing if simulation time was passed
            if (((float)counter) * dt >= simtime)
            {
                if (grid_rank == 0)
                    printProgress(1.);
                break;
            }

            // Writing periodically the progress bar
            if (grid_rank == 0 && counter % 5 == 0)
                printProgress(((float)counter) * dt / simtime);

            // HDF5 writing
            writeChargeHDF5(charge, counter);
            writeExHDF5(Ex, counter);
            writeEyHDF5(Ey, counter);
            writePartHDF5(spec_vec, counter);

            setTime(hdf5_time);

            int flags_coords_mpi[5] = {grid_rank, grid_top,
                                       grid_bottom, grid_right,
                                       grid_left};

            charge->setValue(0.f);

            for (int i = 0; i < Nspecies; i++)
            {
                particle_pusher(Ex, Ey, &spec_vec[i]); // Push particles

                // Loop for sending out of bounds particles (if needed, even to non
                // adjacent processes)
                while (spec_vec[i].advance_cell(flags_coords_mpi))
                {

                    spec_vec[i].prepare_buffer();
                    exchange_particles_buffers(&(spec_vec[i]));
                    spec_vec[i].update_part_list();
                }

                setTime(particle_time);

                get_charge(charge, &spec_vec[i]);

                setTime(field_time);
            }

            exchange_charge_buffers(charge);
            jacobi(phi, charge);
            set_E_value(phi, Ex, Ey);

            setTime(field_time);
        }

        printTime(filename); // Print time info

        closeHDF5(filename);
    }
}