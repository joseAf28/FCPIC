#ifndef _SIMULATION_
#define _SIMULATION_

#include "mpi.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include "field.h"
#include "domain.h"

namespace simulation
{
// Directions
#define RIGHT i + 1
#define LEFT i - 1
#define UP i + Nx
#define DOWN i - Nx
#define P i

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

    // MPI variables - broaden the variables scope
    extern int rank;                                         // rank of the current process
    extern int grid_rank;                                    // rank of the current proces in the virtual grid
    extern int grid_top, grid_bottom, grid_left, grid_right; // ranks of the neighbouring processes
    extern int n_Procs;                                      // total number of processes
    extern int grid[2];                                      // virtual grid dimensions
    extern MPI_Datatype exchange_buffer_type[2];             // MPI_datatype for exchange of buffer cell data
    extern MPI_Comm grid_comm;                               // grid COMMUNICATOR
    extern int offset[2];                                    // offset for cell numbering for subdomains
    extern int grid_coord[2];                                // coordinates of the process in the virtual grid
    extern MPI_Status status;

    // Creates a virtual cartesian topology
    void setup_proc_grid();

    // Packs data for communication across processes
    void setup_MPI_datatypes(int, int, int);

    // sets flags for ghost and buffer cells
    void set_ghost_buffer_flag(domain &);

    // Sets boundary condition values to boundary cells
    void set_bc(field *);

    // exchanges data between processes
    void exchange_buffers(field *, int, int);

    // Jacobi solver
    void jacobi(field *, int, int, field *);

    void write_output_u(domain &, int, int);
    void write_output_charge(domain &, int, int);

}
#endif