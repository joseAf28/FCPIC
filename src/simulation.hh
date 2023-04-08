#ifndef _SIMULATION_
#define _SIMULATION_

#include "mpi.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include "species.hh"
#include "field.hh"
#include "FCPIC_base.hh"

namespace FCPIC
{
    class simulation : public FCPIC_base
    {
    public:
        simulation(int, char **);
        ~simulation() override;

        void readArgs(int, char **);
        void getParamsfromFile(std::string, std::vector<bool>*);
        void setParams();
        void printHelp();
        std::string print_SI(double);
        void printTitle();
        void printProgress(float);

        // Creates a virtual cartesian topology and creates MPI Datatypes
        void setup_proc_grid();
        void get_diagonal_rank(int *, int &);

        // exchanges data between processes
        void exchange_phi_buffers(field *);
        void exchange_charge_buffers(field *);
        void exchange_particles_buffers(species *);

        // Jacobi solver
        void get_charge(field *, species *);
        void jacobi(field *, field *);

        void set_E_value(field *, field *, field *);

        void field_interpolate(field *, field *, float &, float &, part *);
        void init_pusher(field *, field *, species *);
        void particle_pusher(field *, field *, species *);

        // MPI variables
        int grid_rank, rank;                              // rank of the current proces in the virtual grid
        int grid_top, grid_bottom, grid_left, grid_right; // ranks of the neighbouring processes
        int grid_ne, grid_se, grid_nw, grid_sw;           // ranks of diagonal processes: NE, SE, NW, SW

    private:
        MPI_Datatype exchange_field_type[2]; // MPI_datatype for exchange of buffer cell data
        MPI_Comm grid_comm;                  // grid COMMUNICATOR
        int offset[2];                       // offset for cell numbering for subdomains
        int wrap_around[2];
        MPI_Status status;

        // MPI_Datatype exchange_part_type;
        MPI_Aint offsets[8]; // it evaluates to the offset (in bytes) of a given member within a struct or union type
        const int nitems = 8;
        MPI_Datatype exchange_part_type;

        // Simulation variables
        int Nspecies;
        std::vector<int> Npart;
        std::vector<double> charge, mass, temp, vxfluid, vyfluid;
        double aspect, xlen;
        double *X_guard_data, *Y_guard_data;
    };
}
#endif