#ifndef __species__
#define __species__

//FCPIC - 2D Particle-in-Cell code using MPI
//Guilherme Crispim, João Palma, José Afonso
//Advanced Topics in Computational Physics, 2023, IST

//File species.hh:
//Declaration of class Species

#include <string>
#include "field.hh"
#include "FCPIC_base.hh"

////////////////////////////////////////////////////////
// Axis Defintion
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
namespace FCPIC
{
    typedef enum //Flag to check which particles are going be deleted after MPI communication
    {
        BULK,
        SEND
    } Info_part;

    //Struct Particle
    //Represents all physical info of a single particle:
    //-> ix and iy are the indexes of the bottom left grid point nearest to the particle~
    //-> x and y are the particle coordinates in reference to the point of index (ix,iy)
    //(hence x and y are between [0,dx[ and [0,dy[, respectively)
    //-> ux and uy are the particle velocity components
    //-> flag is a flag testing if a particle is part of the process or if must be sent
    typedef struct Particle
    {
        int ix; // Particle cell index
        int iy; // Particle cell index

        float x; // x position inside cell
        float y; // y position inside cell

        float ux;
        float uy;

        bool flag;

    } part;

    //Class species
    //Includes all memory info of the particles inside the process domain.
    //Its constructors are responsible for initializing the particles according to
    //the provided params. It also includes auxiliary functions for preparing comms.
    //It derives from a provided FCPIC_base object to make the simulation params
    //consistent 
    class species : public FCPIC_base
    {
    public:
        species(float, float, float, float *, int *, FCPIC_base const *);
        species(float, float, float, float *, int, FCPIC_base const *);
        ~species() override;

        //Methods used for MPI communication
        void prepare_buffer();
        void update_part_list();
        int advance_cell(int *);

        //Array of particles
        std::vector<part> vec;

        //Buffers to send data from MPI's data exchange
        std::vector<part> send_buffer_north;
        std::vector<part> send_buffer_south;
        std::vector<part> send_buffer_east;
        std::vector<part> send_buffer_west;

        std::vector<part> send_buffer_ne;
        std::vector<part> send_buffer_se;
        std::vector<part> send_buffer_nw;
        std::vector<part> send_buffer_sw;

        //Buffers to receive data from MPI's data exchange
        std::vector<part> recv_buffer_north;
        std::vector<part> recv_buffer_south;
        std::vector<part> recv_buffer_east;
        std::vector<part> recv_buffer_west;

        std::vector<part> recv_buffer_ne;
        std::vector<part> recv_buffer_se;
        std::vector<part> recv_buffer_nw;
        std::vector<part> recv_buffer_sw;

        //Variables to define the size of the arrays for the MPI communication
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

        int np, np_sim; //Total number of particles in the process and in the simulation

        const float q, m;
    };
}
#endif