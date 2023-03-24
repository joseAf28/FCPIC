#include "simulation.hh"
#include "species.hh"

int main(int argc, char **argv)
{
    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::normal_distribution<double> norm(3, 0.5);

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

    // differentiate vector of  each rank
    if (sim.grid_rank == 0) // 0
    {
        vf[0] = 0.;
        vf[1] = 1.7;
        vf[2] = 0;
    }
    if (sim.grid_rank == 1) // 1
    {
        vf[0] = 0;
        vf[1] = 1.6;
        vf[2] = 0.;
    }
    if (sim.grid_rank == 2) // 2
    {
        vf[0] = 0;
        vf[1] = 1.7;
        vf[2] = 0.;
    }
    if (sim.grid_rank == 3) // 3
    {
        vf[0] = 0.;
        vf[1] = 1.8;
        vf[2] = 0.;
    }

    //! constant field for now
    float Ex = 0.;
    float Ey = 0.;

    int flags_coords_mpi[5] = {sim.grid_rank, sim.grid_top, sim.grid_bottom, sim.grid_right, sim.grid_left};

    species test(name, ppc, range, vf, vth);
    test.set_x();
    test.set_u();

    // first iteration of the whole domain
    test.init_pusher(Ex, Ey);
    test.particle_pusher(Ex, Ey);
    test.advance_cell(flags_coords_mpi);

    //
    test.write_output_vec(3.5, sim.grid_rank);
    test.prepare_to_buffer();

    sim.exchange_particles_buffers(&test);

    test.write_input_buffer(4, sim.grid_rank);
    test.write_output_buffer(4, sim.grid_rank);
    test.update_part_list();
    test.write_output_vec(4, sim.grid_rank);

    return 0;
}
