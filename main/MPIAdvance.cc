#include "simulation.hh"
#include "species.hh"
#include <unistd.h>
#include "hdf5.h"

int main(int argc, char **argv)
{
    // initializaing simulation fields and MPI
    FCPIC::simulation *sim = new FCPIC::simulation(argc, argv);
    sim->set_conductive_field_bc();

    // declaring species object
    std::string name = "electron";
    int ppc[2] = {1, 1};
    int range[2] = {20, 10}; // number of cells in each direction

    float *vfa = new float[3];
    float *vfb = new float[3];
    float vth[3] = {0., 0., 0.};
    vfa[0] = 0.;
    vfa[1] = 0.3;
    vfa[2] = 0.;
    vfb[0] = 0.;
    vfb[1] = 0.;
    vfb[2] = 0.;

    // // differentiate vectors
    if (sim->grid_rank == 0) // 0
    {
        vfa[0] = 0.;
        vfa[1] = 0.5;
        vfa[2] = 0;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0;
    }
    if (sim->grid_rank == 1) // 1
    {
        vfa[0] = 0.;
        vfa[1] = 0.5;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }
    if (sim->grid_rank == 2) // 2
    {
        vfa[0] = 0.;
        vfa[1] = 0.3;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }
    if (sim->grid_rank == 3) // 3
    {
        vfa[0] = 0.;
        vfa[1] = 0.3;
        vfa[2] = 0.;
        vfb[0] = 0.;
        vfb[1] = 0.;
        vfb[2] = 0.;
    }

    // fields definition
    FCPIC::field *Ex = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *Ey = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *charge = new FCPIC::field(range[0] + 1, range[1] + 1); // intialize to zero in all entries
    FCPIC::field *phi = new FCPIC::field(range[0] + 1, range[1] + 1);

    // initializing species
    int nb_spec = 2;
    float q[2] = {1, -0.9};

    std::vector<species> spec_vec;

    for (int i = 0; i < nb_spec; i++)
    {
        species test(name, ppc, range, vfa, vth, q[i]);
        spec_vec.push_back(test);
    }

    for (int i = 0; i < nb_spec; i++)
    {
        spec_vec[i].set_x();
        spec_vec[i].set_u();
        spec_vec[i].get_charge(charge);
    }

    sim->jacobi(phi, charge);
    sim->set_E_value(phi, Ex, Ey);

    //! HDF5 Initialization
    std::string h5_name = "../results/test_rank_" + std::to_string(sim->grid_rank) + ".h5";
    const char *h5_char = h5_name.c_str();

    hid_t file_id, group_id, dataset_id, dataspace_id, group; /* identifiers */
    hsize_t dims[2];
    herr_t status;

    file_id = H5Fcreate(h5_char, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = phi->N_x;
    dims[1] = phi->N_y;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    //!!!

    // sim->exchange_charge_buffers(charge);
    // first iteration of the particle pusher
    for (int i = 0; i < nb_spec; i++)
        spec_vec[i].init_pusher(Ex, Ey);

    // // Old Writting
    // std::fstream charge_file;
    // std::string charge_filename = "../results/charge_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(0) + ".txt";
    // charge_file.open(charge_filename, std::ios::out);
    // charge->print_field(charge_file);
    // charge_file.close();
    // //

    for (int counter = 0; counter < 700; counter++)
    {
        // ! Writting in H5 file;
        std::string Ey_name = "Ey_count_" + std::to_string(counter);
        const char *Ey_char = Ey_name.c_str();
        dataset_id = H5Dcreate2(file_id, Ey_char, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(Ey->val[0]));

        std::string Ex_name = "Ex_count_" + std::to_string(counter);
        const char *Ex_char = Ex_name.c_str();
        dataset_id = H5Dcreate2(file_id, Ex_char, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(Ex->val[0]));

        std::string charge_name = "charge_count_" + std::to_string(counter);
        const char *charge_char = charge_name.c_str();
        dataset_id = H5Dcreate2(file_id, charge_char, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(charge->val[0]));
        //!

        // // Writting in file;
        // std::fstream Ex_file;
        // std::string Ex_filename = "../results/Ex_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(counter) + ".txt";
        // Ex_file.open(Ex_filename, std::ios::out);
        // Ex->print_field(Ex_file);
        // Ex_file.close();

        // std::fstream Ey_file;
        // std::string Ey_filename = "../results/Ey_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(counter) + ".txt";
        // Ey_file.open(Ey_filename, std::ios::out);
        // Ey->print_field(Ey_file);
        // Ey_file.close();

        // std::fstream charge_file;
        // std::string charge_filename = "../results/charge_field/rank:_" + std::to_string(sim->grid_rank) + "_counter_" + std::to_string(counter + 1) + ".txt";
        // charge_file.open(charge_filename, std::ios::out);
        // charge->print_field(charge_file);
        // charge_file.close();

        // for (int i = 0; i < nb_spec; i++)
        //     spec_vec[i].write_output_vec(sim->grid_rank, i, counter);
        // //!

        int flags_coords_mpi[5] = {sim->grid_rank, sim->grid_top, sim->grid_bottom, sim->grid_right, sim->grid_left};

        charge->setValue(0.f);

        for (int i = 0; i < nb_spec; i++)
        {
            spec_vec[i].particle_pusher(Ex, Ey);
            spec_vec[i].advance_cell(flags_coords_mpi);
            spec_vec[i].prepare_buffer();

            sim->exchange_particles_buffers(&(spec_vec[i]));

            spec_vec[i].update_part_list();
            spec_vec[i].get_charge(charge);
        }

        // //!jacobi with all the species charge
        sim->jacobi(phi, charge);
        sim->set_E_value(phi, Ex, Ey);
    }
    // std::cout << "End Loop" << std::endl;

    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

    delete Ex;
    delete Ey;
    delete charge, phi;
    delete vfa;
    delete vfb;
    spec_vec.clear();
    delete sim;

    return 0;
}
