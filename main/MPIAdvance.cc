#include "simulation.hh"
#include "species.hh"
#include <unistd.h>

int main(int argc, char **argv)
{
    // initializaing simulation fields and MPI
    FCPIC::simulation *sim = new FCPIC::simulation(argc, argv);
    
    // declaring species object
    std::string name = "electron";
    int ppc[2] = {1, 1};
    int range[2] = {20, 10}; // number of cells in each direction

    float *vfa = new float[3];
    float *vfb = new float[3];
    float temp = 0.;
    vfa[0] = 0.;
    vfa[1] = 0.9;
    vfa[2] = 0.;
    vfb[0] = 0.;
    vfb[1] = 0.;
    vfb[2] = 0.;

    // // differentiate vectors
    if (sim->grid_rank == 0) // 0
    {
        vfa[0] = 0.;
        vfa[1] = 0.19;
        vfa[2] = 0;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0;
    }
    if (sim->grid_rank == 1) // 1
    {
        vfa[0] = 0.;
        vfa[1] = 0.12;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }
    if (sim->grid_rank == 2) // 2
    {
        vfa[0] = 0.;
        vfa[1] = -0.23;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }
    if (sim->grid_rank == 3) // 3
    {
        vfa[0] = 0.;
        vfa[1] = -0.3;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }

    // fields definition
    FCPIC::field *Ex = new FCPIC::field(sim);
    FCPIC::field *Ey = new FCPIC::field(sim);
    FCPIC::field *charge = new FCPIC::field(sim); // intialize to zero in all entries
    FCPIC::field *phi = new FCPIC::field(sim);

    // initializing species
    int nb_spec = 2;
    float q[2] = {1, -0.9};

    std::vector<FCPIC::species> spec_vec;

    for (int i = 0; i < nb_spec; i++)
    {
        FCPIC::species test(name, q[i], 1., temp, vfa, ppc, sim);
        spec_vec.push_back(test);
    }

    for (int i = 0; i < nb_spec; i++)
    {
        sim->get_charge(charge, &spec_vec[i]);
    }

    sim->exchange_charge_buffers(charge);
    

    std::fstream charge_file;
    std::string charge_filename = "../results/charge_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(0) + ".txt";
    charge_file.open(charge_filename, std::ios::out);
    charge->print_field(charge_file);
    charge_file.close();

    // sim->exchange_charge_buffers(charge);

    sim->jacobi(phi, charge);
    sim->set_E_value(phi, Ex, Ey);

    // first iteration of the particle pusher
    for (int i = 0; i < nb_spec; i++)
        sim->init_pusher(Ex, Ey, &spec_vec[i]);
    //spec_vec[i].init_pusher(Ex, Ey);

    for (int counter = 0; counter < 700; counter++)
    {
        // std::cout << "grid_rank: " << sim->grid_rank << "counter:" << counter << std::endl;
        if(sim->grid_rank == 0)
            sim->printProgress(((float) counter)/700.);
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

        for (int i = 0; i < nb_spec; i++)
            spec_vec[i].write_output_vec(sim->grid_rank, i, counter);
        // //!

        int flags_coords_mpi[5] = {sim->grid_rank, sim->grid_top, sim->grid_bottom, sim->grid_right, sim->grid_left};

        charge->setValue(0.f);

        for (int i = 0; i < nb_spec; i++)
        {
            sim->particle_pusher(Ex, Ey, &spec_vec[i]);
            /*
            spec_vec[i].advance_cell(flags_coords_mpi);
            spec_vec[i].prepare_buffer();

            sim->exchange_particles_buffers(&(spec_vec[i]));

            spec_vec[i].update_part_list();
            spec_vec[i].get_charge(charge);
            */
            while(spec_vec[i].advance_cell(flags_coords_mpi)){

            spec_vec[i].prepare_buffer();
            sim->exchange_particles_buffers(&(spec_vec[i]));
            spec_vec[i].update_part_list();

            }
            
            sim->get_charge(charge, &spec_vec[i]);
        }
        sim->exchange_charge_buffers(charge);
        // //!jacobi with all the species charge
        sim->jacobi(phi, charge);
        sim->set_E_value(phi, Ex, Ey);
    }
    // std::cout << "End Loop" << std::endl;
    
    delete Ex;
    delete Ey;
    delete charge, phi;
    delete vfa;
    delete vfb;
    spec_vec.clear();
    
    delete sim;
    return 0;
}
