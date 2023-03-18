#ifndef _SIMULATION_
#define _SIMULATION_

#include "mpi.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include "field.hh"
#include "domain.hh"

namespace FCPIC
{
    // Patch
    typedef enum
    {
        INSIDE,
        XMIN,
        XMAX,
        YMIN,
        YMAX
    } PATCH_type;
    extern PATCH_type patch;

class simulation{
    public:
    simulation(int, char**);
    ~simulation();

    // Creates a virtual cartesian topology
    void setup_proc_grid();

    void set_Xperiodic_field_bc();
    void set_Yperiodic_field_bc();
    void set_Xconductive_field_bc();
    void set_Yconductive_field_bc();

    // Packs data for communication across processes
    void setup_MPI_datatypes(int, int, int);

    // sets flags for ghost and buffer cells
    void set_ghost_buffer_flag(domain &);

    // Sets boundary condition values to boundary cells
    void set_bc(field *);

    // exchanges data between processes
    void exchange_phi_buffers(field *);

    // Jacobi solver
    void jacobi(field *, field *);

    void set_E_value(field *, field *, field *);

    void write_output_u(domain &, int, int);
    void write_output_charge(domain &, int, int);

    private:

    // MPI variables
    //int rank;                                         // rank of the current process
    int grid_rank;                                    // rank of the current proces in the virtual grid
    int grid_top, grid_bottom, grid_left, grid_right; // ranks of the neighbouring processes
    int n_Procs;                                      // total number of processes
    int grid[2];                                      // virtual grid dimensions
    MPI_Datatype exchange_buffer_type[2];             // MPI_datatype for exchange of buffer cell data
    MPI_Comm grid_comm;                               // grid COMMUNICATOR
    int offset[2];                                    // offset for cell numbering for subdomains
    int grid_coord[2];                                // coordinates of the process in the virtual grid
    int wrap_around[2];
    MPI_Status status;

    //Simulation variables
    int Npart;
    double aspect, dx, dy;
    int N_int_x, N_int_y;
    int N_x, N_y, N;
    double *north_recv, *south_recv, *west_recv, *east_recv;
    double *north_send, *south_send, *west_send, *east_send;
    int bc[2];

};
}
#endif