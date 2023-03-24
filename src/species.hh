#ifndef __species__
#define __species__

#include <vector>
#include <string>
#include <memory>
#include <random>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "field.hh"

////////////////////////////////////////////////////////
// Axis Defintion (same conention inside the cell and in the domain as well)
//  y
//  ^
//  |
//  |   .
//  |   .
//  |   .                   ...
//  | (1,0) | (1,1) | (1,2) ...
//  | (0,0) | (0,1) | (0,2) ...
//  0 ------------------------------> x
////////////////////////////////////////////////////////

typedef enum // flag to check which particles are going be deleted after MPI communication
{
    BULK,
    SEND
} Info_part;

//! NO guard cell in the domain of the particles. We only need guard cells in the domain of the charge fields
typedef struct Particle
{
    int ix; // Particle cell index
    int iy; // Particle cell index

    float x; // x position inside cell
    float y; // y position inside cell

    float ux;
    float uy;
    float uz;

    bool flag;

} part;

class species
{
public:
    species(std::string, int *, int *, float *, float *);
    ~species();

    void set_x();
    void set_u();

    void get_charge();

    // methods used for MPI communication
    // !Integrate these methods with the particle pusher
    void prepare_to_buffer();
    void update_part_list();
    void advance_cell(int *);

    // particle pusher - leap frog method
    void init_pusher(const float, const float);
    void particle_pusher(const float, const float);

    //!!!!!!!!!!!!!!!! methods for debugging
    void print();
    void write_output_vec(const int, const int);
    void write_output_buffer(const int, const int);
    void write_input_buffer(const int, const int);

    //!!! dummy
    void update_part();
    void to_buffer();

    // array of particles
    std::vector<part> vec;

    // charge field: Used for the jacobi iteration
    FCPIC::field *charge;

    // buffers to send data from MPI's data exchange
    std::vector<part> send_buffer_north;
    std::vector<part> send_buffer_south;
    std::vector<part> send_buffer_east;
    std::vector<part> send_buffer_west;

    std::vector<part> send_buffer_ne;
    std::vector<part> send_buffer_se;
    std::vector<part> send_buffer_nw;
    std::vector<part> send_buffer_sw;

    // buffers to receive data from MPI's data exchange
    std::vector<part> recv_buffer_north;
    std::vector<part> recv_buffer_south;
    std::vector<part> recv_buffer_east;
    std::vector<part> recv_buffer_west;

    std::vector<part> recv_buffer_ne;
    std::vector<part> recv_buffer_se;
    std::vector<part> recv_buffer_nw;
    std::vector<part> recv_buffer_sw;

    int np; // total number of particles in the simulation

    // simulation box info
    int N_x; // number of x grid points
    int N_y; // number of y grid points

private:
    // species name
    std::string name;

    // mass and charge
    float m = 1; //!! to define in the constructor later
    float q = 1;

    // number of particles per cell
    int ppc[2];

    // number of cells in each direction
    int range[2]; // range[0] -> nb cells in x direction
                  // range[1] -> nb cells in y direction

    // initial general particles momentum
    float vf[3];  // initial fluid velocity
    float vth[3]; // inital thermal velocity

    //!! dx and dy set to 1.0 for now
    float dx = 1.0; // x grid cell size
    float dy = 1.0; // y grid cell size

    float dt = 1.; // time-step - change later

    // Generator of random numbers the thermal boltzmann distribution
    std::mt19937_64 rng;
    std::normal_distribution<double> rand_gauss;
};

#endif