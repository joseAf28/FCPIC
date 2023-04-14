// FCPIC - 2D Particle-in-Cell code using MPI
// Guilherme Crispim, João Palma, José Afonso
// Advanced Topics in Computational Physics, 2023, IST

// File simulation_hdf5.cc:
// Implementation of all HDF5 related functions in the class
// Simulation

#include "simulation.hh"
#include <fstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <unistd.h>

namespace FCPIC
{
    /*
    Function setupHDF5
    Function setupHDF5
    Inputs: source name of the hdf5 output file created. For each rank is created as an output file
    (the rank number is printed in the name) that stores the data from its domain.
    + The hdf5 structure of each file is as follows:
    - group charge where the charge field is recorded at each timestep;
    - group Ex and the group Ey where the x and y  field components are stored at each iteration;
    - group  part_i (i is the species counter) where the struct Part information of each particle of the species i is stored
    - group rank where the grid coordinates of the  MPI virtual topology are stored
    */

    void simulation::setupHDF5(std::string filename)
    {
        std::string h5_name = "../results/" + filename + "_rank_" + std::to_string(grid_rank) + ".h5";
        const char *h5_char = h5_name.c_str();

        // datatype creation to handle the struct Part
        part_id = H5Tcreate(H5T_COMPOUND, sizeof(FCPIC::part));
        H5Tinsert(part_id, "ix", HOFFSET(FCPIC::part, ix), H5T_NATIVE_INT);
        H5Tinsert(part_id, "iy", HOFFSET(FCPIC::part, iy), H5T_NATIVE_INT);
        H5Tinsert(part_id, "x", HOFFSET(FCPIC::part, x), H5T_NATIVE_FLOAT);
        H5Tinsert(part_id, "y", HOFFSET(FCPIC::part, y), H5T_NATIVE_FLOAT);
        H5Tinsert(part_id, "ux", HOFFSET(FCPIC::part, ux), H5T_NATIVE_FLOAT);
        H5Tinsert(part_id, "uy", HOFFSET(FCPIC::part, uy), H5T_NATIVE_FLOAT);
        H5Tinsert(part_id, "flag", HOFFSET(FCPIC::part, flag), H5T_NATIVE_INT);

        file_field = H5Fcreate(h5_char, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t dims[1];
        dims[0] = N_x * N_y;

        if (grid_rank == 0)
        {
            std::cout << "N_x: " << N_x << " N_y: " << N_y << std::endl;
        }

        dataspace_field = H5Screate_simple(1, dims, nullptr);

        group_creation_plist = H5Pcreate(H5P_GROUP_CREATE);
        status_h5 = H5Pset_link_creation_order(group_creation_plist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

        group_charge = H5Gcreate(file_field, "/charge", H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);
        group_Ex = H5Gcreate(file_field, "/Ex", H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);
        group_Ey = H5Gcreate(file_field, "/Ey", H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);

        for (int i = 0; i < Nspecies; i++)
        {
            std::string h5_vec_name = "/part_" + std::to_string(i);
            const char *h5_vec_char = h5_vec_name.c_str();
            hid_t group_idaux = H5Gcreate(file_field, h5_vec_char, H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);
            h5_vec_group.push_back(group_idaux);
        }

        // write the cartesian rank id in the hdf5 file
        hsize_t size_id_rank = 2;
        dataspace_rank = H5Screate_simple(1, &size_id_rank, nullptr);
        group_rank = H5Gcreate(file_field, "/rank", H5P_DEFAULT, group_creation_plist, H5P_DEFAULT);
        std::string rank_id_name = "rank_id";
        const char *rank_id_char = rank_id_name.c_str();
        dataset_rank = H5Dcreate2(group_rank, rank_id_char, H5T_NATIVE_INT, dataspace_rank, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status_h5 = H5Dwrite(dataset_rank, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(grid_coord[0]));

        status_h5 = H5Sclose(dataspace_rank);
        status_h5 = H5Gclose(group_rank);
        status_h5 = H5Dclose(dataset_rank);
    }

    /*
    Function closeHDF5
    + Closes all the HDF5 objects created
    */

    void simulation::closeHDF5(std::string filename)
    {
        hid_t status;

        status = H5Gclose(group_charge);
        status = H5Gclose(group_Ex);
        status = H5Gclose(group_Ey);

        for (int i = 0; i < Nspecies; i++)
            status = H5Gclose(h5_vec_group[i]);

        status = H5Tclose(part_id);
        status = H5Sclose(dataspace_field);
        status = H5Dclose(dataset_field);
        status = H5Fclose(file_field);

        status = H5Pclose(group_creation_plist);

        h5_vec_group.clear();

        if (grid_rank == 0)
        {
            std::cout << "Data written to HDF5 files " << filename << "_rank_[N].h5\n"
                      << std::endl;
        }
    }

    /*
    Functions writeChargeHDF5, writeExHDF5, writeEyHDF5, writePartHDF5
    Inputs: pointer to the object that stores the relevant quantity to store and the number of the simulation's iteration
    + The number of the iteration is used to create distinct datasets that are easy to distinguish
    */
    void simulation::writeChargeHDF5(field *charge, int counter)
    {
        std::string charge_name = "charge_count_" + std::to_string(counter);
        const char *charge_char = charge_name.c_str();
        dataset_field = H5Dcreate2(group_charge, charge_char, H5T_NATIVE_FLOAT, dataspace_field, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status_h5 = H5Dwrite(dataset_field, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(charge->val[0]));
    }
    void simulation::writeExHDF5(field *Ex_field, int counter)
    {
        std::string Ex_name = "Ex_count_" + std::to_string(counter);
        const char *Ex_char = Ex_name.c_str();
        dataset_field = H5Dcreate2(group_Ex, Ex_char, H5T_NATIVE_FLOAT, dataspace_field, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status_h5 = H5Dwrite(dataset_field, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(Ex_field->val[0]));
    }

    void simulation::writeEyHDF5(field *Ey_field, int counter)
    {
        std::string Ey_name = "Ey_count_" + std::to_string(counter);
        const char *Ey_char = Ey_name.c_str();
        dataset_field = H5Dcreate2(group_Ey, Ey_char, H5T_NATIVE_FLOAT, dataspace_field, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status_h5 = H5Dwrite(dataset_field, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(Ey_field->val[0]));
    }

    void simulation::writePartHDF5(std::vector<species> &spec_vec, int counter)
    {
        for (int i = 0; i < Nspecies; i++)
        {
            std::string part_name = "part_count_" + std::to_string(counter);
            const char *part_char = part_name.c_str();
            hsize_t vec_size = spec_vec[i].vec.size();
            dataspace_part = H5Screate_simple(1, &vec_size, nullptr);
            dataset_part = H5Dcreate2(h5_vec_group[i], part_char, part_id, dataspace_part, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status_h5 = H5Dwrite(dataset_part, part_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(spec_vec[i].vec[0]));

            H5Dclose(dataset_part);
            H5Sclose(dataspace_part);
        }
    }
}