#include "hdf5.h"
#include "simulation.hh"

int main()
{

    // hid_t file_id, group_id, dataset_id, dataspace_id, group; /* identifiers */
    // hsize_t dims[2];
    // herr_t status;

    // std::vector<double> vec1(100, 7.);
    // std::vector<double> vec2(100, 8.);

    // int rank = 0;
    // std::string h5_name = "Xdata_rank_" + std::to_string(rank) + ".h5";
    // const char *h5_char = h5_name.c_str();

    // file_id = H5Fcreate(h5_char, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // dims[0] = 2;
    // dims[1] = 10;
    // dataspace_id = H5Screate_simple(2, dims, NULL);

    // dataset_id = H5Dcreate2(file_id, "Ey", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vec1[0]));

    // dataset_id = H5Dcreate2(file_id, "Ex", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vec2[0]));

    // status = H5Sclose(dataspace_id);
    // status = H5Dclose(dataset_id);
    // status = H5Fclose(file_id);

    // write structure hdf5

    typedef struct s1_t
    {
        int a;
        float b;
        double c;
    } s1_t;

    hid_t s1_tid;

    std::vector<s1_t> vec;

    s1_t A_dummy;
    A_dummy.a = 1;
    A_dummy.b = 18.;
    A_dummy.c = 3 / 7;
    vec.assign(10, A_dummy);

    int i;
    hid_t file, dataset, space, group;
    herr_t status;
    hsize_t dim = 4;

    file = H5Fcreate("h5Struct.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    group = H5Gcreate(file, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(s1_t));
    H5Tinsert(s1_tid, "a_name", HOFFSET(s1_t, a), H5T_NATIVE_INT);
    H5Tinsert(s1_tid, "b_name", HOFFSET(s1_t, b), H5T_NATIVE_FLOAT);
    H5Tinsert(s1_tid, "c_name", HOFFSET(s1_t, c), H5T_NATIVE_DOUBLE);

    for (int i = 0; i < 3; i++)
    {
        std::string name = "name_" + std::to_string(i);
        const char *name_char = name.c_str();
        space = H5Screate_simple(1, &dim, NULL);
        dataset = H5Dcreate2(group, name_char, s1_tid, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vec[0]));
        H5Dclose(dataset);
        H5Sclose(space);
        dim = dim + 1;
    }

    H5Tclose(s1_tid);
    H5Gclose(group);
    // H5Sclose(space);
    // H5Dclose(dataset);
    H5Fclose(file);

    return 0;
}