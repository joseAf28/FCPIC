#include "simulation.hh"
#include "species.hh"
#include <unistd.h>

int main(int argc, char **argv)
{
    // not used for now
    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::normal_distribution<double> norm(3, 0.5);

    // initializaing simulation fields and MPI
    FCPIC::simulation *sim = new FCPIC::simulation(argc, argv);
    sim->set_conductive_field_bc();

    // declaring species object
    std::string name = "electron";
    int ppc[2] = {1, 1};
    int range[2] = {10, 10}; // number of cells in each direction

    float *vfa = new float[3];
    float *vfb = new float[3];
    float vth[3] = {0., 0., 0.};
    vfa[0] = 0.9;
    vfa[1] = 0.;
    vfa[2] = 0.;
    vfb[0] = -0.;
    vfb[1] = 0.;
    vfb[2] = 0.;

    // // differentiate vectors
    if (sim->grid_rank == 0) // 0
    {
        vfa[0] = 0.6;
        vfa[1] = 0.;
        vfa[2] = 0;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0;
    }
    if (sim->grid_rank == 1) // 1
    {
        vfa[0] = 0.5;
        vfa[1] = 0.;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }
    if (sim->grid_rank == 2) // 2
    {
        vfa[0] = -0.9;
        vfa[1] = 0.;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }
    if (sim->grid_rank == 3) // 3
    {
        vfa[0] = -0.3;
        vfa[1] = 0.;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }

    FCPIC::field *Ex = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *Ey = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *charge;
    FCPIC::field *phi = new FCPIC::field(range[0] + 1, range[1] + 1);

    // Ex->setValue(1.0);
    // Ey->setValue(1.5);

    // initializing species
    species specA(name, ppc, range, vfa, vth, 1.);
    species specB(name, ppc, range, vfa, vth, -0.9);
    specA.set_x();
    specA.set_u();

    specB.set_x();
    specB.set_u();

    // if (sim->grid_rank == 0)
    // {
    //     for (int i = 1; i < specA.vec.size(); i++)
    //         specA.vec[i].flag = SEND;

    //     specA.vec.erase(std::remove_if(specA.vec.begin(), specA.vec.end(), [&specA](const part obj)
    //                                    { return (obj.flag == SEND); }),
    //                     specA.vec.end());

    //     std::cout << sim->grid_rank << "vec.size: " << specA.vec.size() << std::endl;
    // }

    // if (sim->grid_rank != 0)
    // {
    //     for (int i = 0; i < specA.vec.size(); i++)
    //         specA.vec[i].flag = SEND;

    //     specA.vec.erase(std::remove_if(specA.vec.begin(), specA.vec.end(), [&specA](const part obj)
    //                                    { return (obj.flag == SEND); }),
    //                     specA.vec.end());
    //     std::cout << sim->grid_rank << "vec.size: " << specA.vec.size() << std::endl;
    // }

    // specA.vec[0].ix = 3;
    // specA.vec[0].iy = 3;
    specA.write_output_vec(-1, sim->grid_rank);

    // std::cout << "grid: " << sim->grid_rank << " vec.sizeA: " << specA.vec.size() << " vec.size B: " << specB.vec.size() << std::endl;

    charge = new FCPIC::field(range[0] + 1, range[1] + 1);
    specA.get_charge(charge); // getting initial charge field
    // specB.get_charge(charge);

    std::fstream charge_file;
    std::string charge_filename = "../results/charge_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(0) + ".txt";
    charge_file.open(charge_filename, std::ios::out);
    charge->print_field(charge_file);
    charge_file.close();

    // sim->exchange_charge_buffers(charge);

    sim->jacobi(phi, charge);
    sim->set_E_value(phi, Ex, Ey);

    // getting first consistent electric field
    // if(sim.grid_rank!=0)
    //    sleep(3);
    // sleep(sim.grid_rank);

    // phi->print_field(std::cout);

    specA.init_pusher(Ex, Ey); // first iteration of the particle pusher
    // specB.init_pusher(Ex, Ey);
    // specA.size(sim->grid_rank);
    specA.update_part_list();

    for (int counter = 0; counter < 200; counter++)
    {
        // std::cout << "grid_rank: " << sim->grid_rank << "counter:" << counter << std::endl;
        // //! Writting in file;
        std::fstream Ex_file;
        std::string Ex_filename = "../results/Ex_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(counter) + ".txt";
        Ex_file.open(Ex_filename, std::ios::out);
        Ex->print_field(Ex_file);
        Ex_file.close();

        std::fstream Ey_file;
        std::string Ey_filename = "../results/Ey_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(counter) + ".txt";
        Ey_file.open(Ey_filename, std::ios::out);
        Ey->print_field(Ey_file);
        Ey_file.close();

        std::fstream charge_file;
        std::string charge_filename = "../results/charge_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(counter + 1) + ".txt";
        charge_file.open(charge_filename, std::ios::out);
        charge->print_field(charge_file);
        charge_file.close();
        // specA.write_output_vec(counter, sim->grid_rank); // debugging

        int flags_coords_mpi[5] = {sim->grid_rank, sim->grid_top, sim->grid_bottom, sim->grid_right, sim->grid_left};

        specA.particle_pusher(Ex, Ey); // includes field interpolation at particle position
        // while (specA.advance_cell(flags_coords_mpi))
        // {
        //     // test.write_output_vec(counter, sim.grid_rank); // debugging
        //     // specA.size(sim->grid_rank);
        //     specA.prepare_buffer();
        //     std::cout << "grid: " << sim->grid_rank << " vec.sizeA: " << specA.vec.size() << " vec.size B: " << specB.vec.size() << std::endl;

        //     sim->exchange_particles_buffers(&specA);

        //     // test.write_input_buffer(counter, sim.grid_rank);  // debugging communication

        //     specA.update_part_list(); // update the list of particles of each ptocess because of the MPI exchange
        //     std::cout << "****************" << std::endl;
        // }
        specA.advance_cell(flags_coords_mpi);
        specA.prepare_buffer();

        sim->exchange_particles_buffers(&specA);
        specA.update_part_list();
        // //specB.particle_pusher(Ex, Ey);

        // while (specB.advance_cell(flags_coords_mpi))
        // {
        //     specB.prepare_buffer();
        //     sim->exchange_particles_buffers(&specB);

        //     specB.update_part_list();
        //     std::cout << "****************" << std::endl;
        // }

        // getting charge distribution and doing the sum at the same time
        charge->setValue(0.f);

        specA.get_charge(charge);
        // specB.get_charge(charge);

        // jacobi with all the species charge
        sim->jacobi(phi, charge);
        sim->set_E_value(phi, Ex, Ey);

        // sleep(sim.grid_rank);
        // charge->print_field(std::cout);

        // test.write_output_vec(counter, sim.grid_rank); // debugging
        //  std::cout << "End grid_rank: " << sim.grid_rank << std::endl;
    }
    // std::cout << "End Loop" << std::endl;

    delete Ex;
    delete Ey;
    delete charge, phi;
    delete vfa;
    delete vfb;
    delete sim;
    return 0;
}
