#include "simulation.hh"

int main(int argc, char **argv)
{
    // initializaing simulation fields and MPI
    FCPIC::simulation *sim = new FCPIC::simulation(argc, argv);

    if(sim->sim_true){
        // declaring species object
        std::string name = "electron";
        int ppc[2] = {1, 1};

        float *vfa = new float[3];

        // fields definition
        FCPIC::field *Ex = new FCPIC::field(sim);
        FCPIC::field *Ey = new FCPIC::field(sim);
        FCPIC::field *charge = new FCPIC::field(sim);
        FCPIC::field *phi = new FCPIC::field(sim);
        std::vector<FCPIC::species> spec_vec;

        // initializing species
        int nb_spec = sim->Nspecies;

        int flags_coords_mpi[5] = {sim->grid_rank, sim->grid_top, 
                                   sim->grid_bottom, sim->grid_right, 
                                   sim->grid_left};

        for (int i = 0; i < nb_spec; i++)
        {
            vfa[0] = sim->vxfluid[i];
            vfa[1] = sim->vyfluid[i];
            //FCPIC::species test(name, sim->charge[i], sim->mass[i], sim->temp[i], vfa, ppc, sim);
            FCPIC::species test(name, sim->charge[i], sim->mass[i], sim->temp[i], vfa, sim->Npart[i], sim);
            spec_vec.push_back(test);
        }
        
        for (int i = 0; i < nb_spec; i++){
        while (spec_vec[i].advance_cell(flags_coords_mpi))
            {
                spec_vec[i].prepare_buffer();
                sim->exchange_particles_buffers(&(spec_vec[i]));
                spec_vec[i].update_part_list();
            }}
        
        sim->run_simulation(Ex, Ey, phi, charge, spec_vec, "final_sim_small");
        

        delete Ex;
        delete Ey;
        delete charge, phi;

        delete vfa;
        spec_vec.clear();
    }

    delete sim;

    return 0;
}
