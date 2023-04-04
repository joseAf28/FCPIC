#include <iostream>
#include <string>
#include "H5Cpp.h"
#include <random>

using namespace H5;
int main()
{

    const H5std_string FILE_NAME("data.h5");
    const H5std_string DATASET_NAME("DOUBLEArray");
    const int NX = 123; // dataset dimensions
    const int NY = 4563;
    const int RANK = 2;

    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::normal_distribution<double> norm(3, 0.5);

    int i, j;
    double data_arr[NX][NY]; // buffer for data to write
    for (j = 0; j < NX; j++)
    {
        for (i = 0; i < NY; i++)
            data_arr[j][i] = norm(rng);
    }

    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    hsize_t dimsf[2]; // dataset dimensions
    dimsf[0] = NX;
    dimsf[1] = NY;
    DataSpace dataspace(RANK, dimsf);
    /*
     * Define datatype for the data in the file.
     * We will store little endian DOUBLE numbers.
     */
    FloatType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);

    dataset.write(data_arr, PredType::NATIVE_DOUBLE);

    return 0;
}

// need to compile with h5c++ instead of g++