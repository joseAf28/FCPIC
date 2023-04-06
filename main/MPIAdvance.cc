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
    int ppc[2] = {2, 2};
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
    // sim->exchange_charge_buffers(charge);

    //! HDF5 Initialization
    std::string h5_name = "../results/data_rank_" + std::to_string(sim->grid_rank) + ".h5";
    const char *h5_char = h5_name.c_str();

    hid_t file_field, dataset_field, dataspace_field;
    hid_t dataspace_part, dataset_part;
    hid_t group_charge, group_Ex, group_Ey, group_particles;
    hsize_t dims[2];
    herr_t status;

    std::vector<hid_t> h5_vec_group;

    hid_t part_id;
    part_id = H5Tcreate(H5T_COMPOUND, sizeof(part));
    H5Tinsert(part_id, "ix", HOFFSET(part, ix), H5T_NATIVE_INT);
    H5Tinsert(part_id, "iy", HOFFSET(part, iy), H5T_NATIVE_INT);
    H5Tinsert(part_id, "x", HOFFSET(part, x), H5T_NATIVE_FLOAT);
    H5Tinsert(part_id, "y", HOFFSET(part, y), H5T_NATIVE_FLOAT);
    H5Tinsert(part_id, "ux", HOFFSET(part, ux), H5T_NATIVE_FLOAT);
    H5Tinsert(part_id, "uy", HOFFSET(part, uy), H5T_NATIVE_FLOAT);
    H5Tinsert(part_id, "uz", HOFFSET(part, uz), H5T_NATIVE_FLOAT);
    H5Tinsert(part_id, "flag", HOFFSET(part, flag), H5T_NATIVE_INT);

    file_field = H5Fcreate(h5_char, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = phi->N_x;
    dims[1] = phi->N_y;
    dataspace_field = H5Screate_simple(2, dims, nullptr);

    hid_t group_creation_plist;
    group_creation_plist = H5Pcreate(H5P_GROUP_CREATE);
    status = H5Pset_link_creation_order(group_creation_plist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

    group_charge = H5Gcreate(file_field, "/charge", H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);
    group_Ex = H5Gcreate(file_field, "/Ex", H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);
    group_Ey = H5Gcreate(file_field, "/Ey", H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);

    for (int i = 0; i < nb_spec; i++)
    {
        std::string h5_vec_name = "/part_" + std::to_string(i);
        const char *h5_vec_char = h5_vec_name.c_str();
        hid_t group_idaux = H5Gcreate(file_field, h5_vec_char, H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);
        h5_vec_group.push_back(group_idaux);
    }
    //!!!

    // first iteration of the particle pusher
    for (int i = 0; i < nb_spec; i++)
        spec_vec[i].init_pusher(Ex, Ey);

    for (int counter = 0; counter < 700; counter++)
    {
        // ! Writting in H5 file;
        std::string Ey_name = "Ey_count_" + std::to_string(counter);
        const char *Ey_char = Ey_name.c_str();
        dataset_field = H5Dcreate2(group_Ey, Ey_char, H5T_NATIVE_DOUBLE, dataspace_field, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(Ey->val[0]));

        std::string Ex_name = "Ex_count_" + std::to_string(counter);
        const char *Ex_char = Ex_name.c_str();
        dataset_field = H5Dcreate2(group_Ex, Ex_char, H5T_NATIVE_DOUBLE, dataspace_field, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(Ex->val[0]));

        std::string charge_name = "charge_count_" + std::to_string(counter);
        const char *charge_char = charge_name.c_str();
        dataset_field = H5Dcreate2(group_charge, charge_char, H5T_NATIVE_DOUBLE, dataspace_field, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(charge->val[0]));

        // loop over all species
        for (int i = 0; i < nb_spec; i++)
        {
            std::string part_name = "part_count_" + std::to_string(counter);
            const char *part_char = part_name.c_str();
            hsize_t vec_size = spec_vec[i].vec.size();
            dataspace_part = H5Screate_simple(1, &vec_size, nullptr);
            dataset_part = H5Dcreate2(h5_vec_group[i], part_char, part_id, dataspace_part, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataset_part, part_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(spec_vec[i].vec[0]));

            H5Dclose(dataset_part);
            H5Sclose(dataspace_part);
        }
        //!

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

    status = H5Gclose(group_charge);
    status = H5Gclose(group_Ex);
    status = H5Gclose(group_Ey);

    for (int i = 0; i < nb_spec; i++)
        status = H5Gclose(h5_vec_group[i]);

    status = H5Tclose(part_id);
    status = H5Sclose(dataspace_field);
    status = H5Dclose(dataset_field);
    status = H5Fclose(file_field);

    status = H5Pclose(group_creation_plist);

    delete Ex;
    delete Ey;
    delete charge, phi;
    delete vfa;
    delete vfb;
    spec_vec.clear();
    h5_vec_group.clear();

    delete sim;

    return 0;
}
