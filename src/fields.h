#ifndef __fields__
#define __fields__

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "species.h"

using namespace std;

class fields
{
public:
    fields(const int, const int, vector<float> &);
    ~fields();

    void potential_solver();                    // uses Jacobi Method to solve the Poisson Equation
    void field_solver();                        // Uses central differences to calculate E field from Potential
    void field_inter(part &, float &, float &); // Interpolates the field inside the cell at particles' position

    void print();

private:
    vector<float> pot_vec;
    vector<float> charge_vec;

    vector<float> Ex_vec;
    vector<float> Ey_vec;

    // simulation box info
    int nx; // number of x grid points
    int ny; // number of y grid points

    //!! dx and dy set to 1.0 for now
    float dx = 1.0; // x grid cell size
    float dy = 1.0; // y grid cell size

    float xbox; // size x axis of simulation box
    float ybox; // size y axis of simulation box

    //! Using Jacobi Iteration Method
    int max_iter = 10000;
    float tol = 1e-6;
};

#endif