#include "simulation.hh"

int main(int argc, char **argv)
{
    // initializaing simulation fields and MPI
    FCPIC::simulation *sim = new FCPIC::simulation(argc, argv);

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

    for (int i = 0; i < nb_spec; i++)
    {
        vfa[0] = sim->vxfluid[i];
        vfa[1] = sim->vyfluid[i];
        FCPIC::species test(name, sim->charge[i], sim->mass[i], sim->temp[i], vfa, ppc, sim);
        spec_vec.push_back(test);
    }

    sim->run_simulation(Ex, Ey, phi, charge, spec_vec);

    hid_t status;

    status = H5Gclose(sim->group_charge);
    status = H5Gclose(sim->group_Ex);
    status = H5Gclose(sim->group_Ey);

    for (int i = 0; i < nb_spec; i++)
        status = H5Gclose(sim->h5_vec_group[i]);

    status = H5Tclose(sim->part_id);
    status = H5Sclose(sim->dataspace_field);
    status = H5Dclose(sim->dataset_field);
    status = H5Fclose(sim->file_field);

    status = H5Pclose(sim->group_creation_plist);

    sim->h5_vec_group.clear();

    delete Ex;
    delete Ey;
    delete charge, phi;

    delete vfa;
    spec_vec.clear();
    delete sim;

    return 0;
}
