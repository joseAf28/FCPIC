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

    int set_nb_part();
    void set_X();
    void set_U();

    void get_current(vector<float> &);
    void advance();

    void print();

private:
    // species name
    string name;

    // array of particles
    unique_ptr<part> vec;

    // mass and charge
    float m = 1; // to define later
    float q = 1;

    // number of particles per cell
    int ppc[2];

    int np; // number of particles

    int range[4];      // range[0] , range[1] -> xmin, xmax-1
                       // range[2] , range[3] -> ymin, ymax-1
    int init_U, end_U; // index where the particles have defined momentum

    // initial general particles momentum
    float vf[3];
    float vth[3];

    // simulation box info
    int nx;    // number of grid points
    float dx;  // grid cell size
    float box; // size of simulation box

    // Generator of random numbers
    mt19937_64 rng;
    normal_distribution<double> rand_gauss;
};

#endif