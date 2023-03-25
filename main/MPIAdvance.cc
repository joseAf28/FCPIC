#include "simulation.hh"
#include "species.hh"

int main(int argc, char **argv)
{
    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::normal_distribution<double> norm(3, 0.5);

    // initializaing simulation fields and MPI
    FCPIC::simulation sim(argc, argv);
    sim.set_conductive_field_bc();

    // declaring species object
    std::string name = "electron";
    int ppc[2] = {1, 1};
    int range[2] = {10, 10}; // number of cells in each direction

    float *vf = new float[3];
    float vth[3] = {0.0, 0.0, 0.0};

    vf[0] = 0.;
    vf[1] = 1.3;
    vf[2] = 0;

    // differentiate vectors
    if (sim.grid_rank == 0) // 0
    {
        vf[0] = 1.8;
        vf[1] = 1.3;
        vf[2] = 0;
    }
    if (sim.grid_rank == 1) // 1
    {
        vf[0] = -1.5;
        vf[1] = -1.6;
        vf[2] = 0.;
    }
    if (sim.grid_rank == 2) // 2
    {
        vf[0] = 1.3;
        vf[1] = -1.1;
        vf[2] = 0.;
    }
    if (sim.grid_rank == 3) // 3
    {
        vf[0] = -1.5;
        vf[1] = -1.7;
        vf[2] = 0.;
    }

    //! constant field for now
    float Ex = 0.;
    float Ey = 0.;

    // initializing species
    species test(name, ppc, range, vf, vth);
    test.set_x();
    test.set_u();

    // getting initial charge field
    // getting first electric field
    test.init_pusher(Ex, Ey); // first iteration of the particle pusher

    for (int counter = 0; counter < 5; counter++)
    {
        int flags_coords_mpi[5] = {sim.grid_rank, sim.grid_top, sim.grid_bottom, sim.grid_right, sim.grid_left};

        // getting charge distribution
        // compute de potential field with the jacobi iteration
        // getting the E field interpolation

        test.particle_pusher(Ex, Ey);
        test.advance_cell(flags_coords_mpi);

        test.write_output_vec(counter, sim.grid_rank); // debugging species list

        test.prepare_buffer();
        // std::cout << "Init grid_rank: " << sim.grid_rank << std::endl;
        sim.exchange_particles_buffers(&test);

        test.write_input_buffer(counter, sim.grid_rank);  // debugging communication
        test.write_output_buffer(counter, sim.grid_rank); // debugging

        test.update_part_list();                       // update the list of particles of each ptocess because of the MPI exchange
        test.write_output_vec(counter, sim.grid_rank); // debugging
        // std::cout << "End grid_rank: " << sim.grid_rank << std::endl;
    }
    std::cout << "End Loop" << std::endl;

    return 0;
}
