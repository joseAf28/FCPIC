#include "hdf5.h"
#include "simulation.hh"
#define FILE "groups.h5"

int main()
{

    hid_t file_id, group_id, dataset_id, dataspace_id, group, datasetGroup_id; /* identifiers */
    hsize_t dims[2];
    herr_t status;
    int i, j, dset1_data[3][3], dset2_data[2][10];

    std::vector<double> vec1(100, 7.);
    std::vector<double> vec2(100, 8.);

    int rank = 0;
    std::string h5_name = "Xdata_rank_" + std::to_string(rank) + ".h5";
    const char *h5_char = h5_name.c_str();
    file_id = H5Fcreate(h5_char, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create a group named "/MyGroup" in the file. */

    // dims[0] = 3;
    // dims[1] = 3;
    // dataspace_id = H5Screate_simple(2, dims, NULL);

    // dataset_id = H5Dcreate2(file_id, "Ex", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vec1[0]));
    // status = H5Sclose(dataspace_id);
    // status = H5Dclose(dataset_id);

    // /* Open an existing group of the specified file. */
    // group_id = H5Gopen2(file_id, "/Data", H5P_DEFAULT);

    /* Create the data space for the second dataset. */
    dims[0] = 2;
    dims[1] = 10;
    // group_id = H5Gcreate(file_id, "MyGroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataspace_id = H5Screate_simple(2, dims, NULL);
    /* Create the second dataset in group "Group_A". */
    // dataset_id = H5Dcreate2(group_id, "dset2", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    datasetGroup_id = H5Dcreate2(file_id, "dset2", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataset_id = H5Dcreate2(file_id, "Ey", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vec2[0]));

    dataset_id = H5Dcreate2(file_id, "Ex", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vec1[0]));

    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);

    // dataset_id = H5Dcreate2(group_id, "dset2", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Close the group. */
    // status = H5Gclose(group_id);

    /* Close the file. */
    status = H5Fclose(file_id);
}