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

    //!! Test: All particles forced to move to top right corner
    //! (external) constant field for now

    FCPIC::field *Ex = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *Ey = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *charge = new FCPIC::field(range[0]+1, range[1]+1);
    FCPIC::field *phi = new FCPIC::field(range[0]+1, range[1]+1);

    //Ex->setValue(1.0);
    //Ey->setValue(1.5);

    // initializing species
    species test(name, ppc, range, vf, vth);
    test.set_x();
    test.set_u();

    test.get_charge(charge); // getting initial charge field
    sim.exchange_charge_buffers(charge);
    sim.jacobi(phi, charge);
    sim.set_E_value(phi, Ex, Ey);
    // getting first consistent electric field
    //if(sim.grid_rank!=0)
    //    sleep(3);
    //sleep(sim.grid_rank);
    
    //phi->print_field(std::cout);

    test.init_pusher(Ex, Ey); // first iteration of the particle pusher
    
    /*
    for (int counter = 0; counter < 50; counter++)
    {
        int flags_coords_mpi[5] = {sim.grid_rank, sim.grid_top, sim.grid_bottom, sim.grid_right, sim.grid_left};
        
        // compute de potential field with the jacobi iteration
        // getting the E field in the grid

        test.particle_pusher(Ex, Ey); // includes field interpolation at particle position
        while(test.advance_cell(flags_coords_mpi)){
            //test.write_output_vec(counter, sim.grid_rank); // debugging

            test.prepare_buffer();
            sim.exchange_particles_buffers(&test);

            test.write_input_buffer(counter, sim.grid_rank);  // debugging communication
            test.write_output_buffer(counter, sim.grid_rank); // debugging

            test.update_part_list(); // update the list of particles of each ptocess because of the MPI exchange
        }

        test.get_charge(charge);       // getting charge distribution
        sim.exchange_charge_buffers(charge);
        sim.jacobi(phi, charge);
        sim.set_E_value(phi, Ex, Ey);
        //sleep(sim.grid_rank);
        //charge->print_field(std::cout);

        //test.write_output_vec(counter, sim.grid_rank); // debugging
        // std::cout << "End grid_rank: " << sim.grid_rank << std::endl;
    }
    std::cout << "End Loop" << std::endl;

    if(sim.grid_rank==2){
        charge->print_field(std::cout);
        std::cout << "\n";
        phi->print_field(std::cout);
        std::cout << "\n";
        Ex->print_field(std::cout);
        std::cout << "\n";
        Ey->print_field(std::cout);
        std::cout << "\n";
    }

    delete Ex;
    delete Ey;
    delete charge, phi;
    delete vf;
    */
    return 0;
}
