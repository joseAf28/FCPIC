#ifndef __species__
#define __species__

#include <vector>
#include <string>
#include <memory>
#include <random>
#include <algorithm>
#include <iostream>

using namespace std;

typedef struct Particle
{
    int ix; // Particle cell index
    int iy; // Particle cell index

    float x; // x position inside cell
    float y; // y position inside cell

    float ux;
    float uy;
    float uz;
} part;

class species
{
public:
    species(string, int *, int *, int *, float *, float *);
    ~species();

    int set_nb();
    void set_X();
    void set_U();

    void get_charge(vector<float> &);
    void advance_cell(int);

    void get_grid_points(int &, int &);

    void print();

    // array of particles
    unique_ptr<part> vec;

private:
    // species name
    string name;

    // mass and charge
    float m = 1; //!! to define in the constructor later
    float q = 1;

    // number of particles per cell
    int ppc[2];

    int np; // number of particles in the simulation

    //(i, j) pairs per cell
    int range[4];      // range[0] , range[1] -> i_min, i_max-1
                       // range[2] , range[3] -> j_min, j_max-1
    int init_U, end_U; // index where the particles have non zero momentum
                       // end_U equals to nb for now

    // initial general particles momentum
    float vf[3];  // initial fluid velocity
    float vth[3]; // inital thermal velocity

    // simulation box info
    int nx; // number of x grid points
    int ny; // number of y grid points

    //!! dx and dy set to 1.0 for now
    float dx = 1.0; // x grid cell size
    float dy = 1.0; // y grid cell size

    float xbox; // size x axis of simulation box
    float ybox; // size y axis of simulation box

    // ?We define the cartesian coordinates as:
    // ?The origin of the cell is at the top left corner of the cell.

    // Generator of random numbers the thermal boltzmann distribution
    mt19937_64 rng;
    normal_distribution<double> rand_gauss;
};

#endif