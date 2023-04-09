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
#include "FCPIC_base.hh"

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
namespace FCPIC{
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

class species : public FCPIC_base
{
public:
    species(std::string, float, float, float, float *, int *, FCPIC_base const *);
    species(std::string, float, float, float, float *, int, FCPIC_base const *);
    ~species() override;

    // methods used for MPI communication
    void prepare_buffer();
    void update_part_list();
    int advance_cell(int *);

    // particle pusher - leap frog method
    //void field_interpolate(field *, field *, float &, float &, part *);
    //void init_pusher(field *, field *);
    //void particle_pusher(field *, field *);

    //!!!!!!!!!!!!!!!! emporary methods for debugging
    void print();
    void write_output_vec(const int, const int, const int);
    void write_output_buffer(const int, const int);
    void write_input_buffer(const int, const int);

    // array of particles
    std::vector<part> vec;

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

    // variables to define the size of the arrays for the MPI communication
    int size_send_north = 0;
    int size_send_south = 0;
    int size_send_east = 0;
    int size_send_west = 0;

    int size_send_ne = 0;
    int size_send_se = 0;
    int size_send_nw = 0;
    int size_send_sw = 0;

    int size_recv_north = 0;
    int size_recv_south = 0;
    int size_recv_east = 0;
    int size_recv_west = 0;

    int size_recv_ne = 0;
    int size_recv_se = 0;
    int size_recv_nw = 0;
    int size_recv_sw = 0;

    int np, np_sim; // total number of particles in the process and in the simulation

    const float q, m;

private:
    // species name
    std::string name;
};
}
#endif