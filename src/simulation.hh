#ifndef _SIMULATION_
#define _SIMULATION_

#include "mpi.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include "species.hh"
#include "field.hh"
#include "domain.hh"

namespace FCPIC
{
    class simulation
    {
    public:
        simulation(int, char **);
        ~simulation();

        void readArgs(int, char **);
        void getParamsfromFile(std::string, std::vector<bool>*);
        void setParams();
        void printHelp();
        std::string print_SI(double);
        void printTitle();

        // Creates a virtual cartesian topology and creates MPI Datatypes
        void setup_proc_grid();
        void get_diagonal_rank(int *, int &);

        // Sets boundary condition values to boundary cells
        void set_periodic_field_bc();
        void set_conductive_field_bc();

        // exchanges data between processes
        void exchange_phi_buffers(field *);
        void exchange_charge_buffers(field *);
        void exchange_particles_buffers(species *);

        // Jacobi solver
        void jacobi(field *, field *);

        void set_E_value(field *, field *, field *);

        // MPI variables
        int grid_rank, rank;                              // rank of the current proces in the virtual grid
        int grid_top, grid_bottom, grid_left, grid_right; // ranks of the neighbouring processes
        int grid_ne, grid_se, grid_nw, grid_sw;           // ranks of diagonal processes: NE, SE, NW, SW
        int wrap_around[2];                               // public for the class species particle pusher can see it

    private:
        int n_Procs;                         // total number of processes
        int grid[2];                         // virtual grid dimensions
        MPI_Datatype exchange_field_type[2]; // MPI_datatype for exchange of buffer cell data
        MPI_Comm grid_comm;                  // grid COMMUNICATOR
        int offset[2];                       // offset for cell numbering for subdomains
        int grid_coord[2];                   // coordinates of the process in the virtual grid
        MPI_Status status;

        // MPI_Datatype exchange_part_type;
        MPI_Aint offsets[8]; // it evaluates to the offset (in bytes) of a given member within a struct or union type
        const int nitems = 8;
        MPI_Datatype exchange_part_type;

        // Simulation variables
        std::vector<int> Npart;
        std::vector<double> charge, mass, temp, vxfluid, vyfluid;
        double Lref, Vref, Tref, Nref, aspect, xlen;
        int Nspecies, nxproc, N_total_x, N_total_y;
        double simtime;
        double dx, dy, dt;
        int N_int_x, N_int_y;
        int N_x, N_y, N;
        double *X_guard_data, *Y_guard_data;
        double *X_guard_data1, *Y_guard_data1;
        double *X_guard_data2, *Y_guard_data2;
        int bc[2];
    };
}
#endif