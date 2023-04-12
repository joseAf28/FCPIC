#ifndef _SIMULATION_
#define _SIMULATION_

#include "mpi.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include "species.hh"
#include "field.hh"
#include "FCPIC_base.hh"
#include "hdf5.h"

namespace FCPIC
{
    class simulation : public FCPIC_base
    {
    public:
        simulation(int, char **);
        ~simulation() override;

        void readArgs(int, char **);
        void getParamsfromFile(std::string, std::vector<bool> *);
        void setParams();
        void confirmParams();
        void printHelp();
        std::string print_SI(double, int);
        void printTitle();
        void printProgress(float);
        void printTime(std::string);

        void setTime(float &);
        void setTime();

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

        void setupHDF5(std::string);
        void closeHDF5();
        void writeChargeHDF5(field *, int);
        void writeExHDF5(field *, int);
        void writeEyHDF5(field *, int);
        void writePartHDF5(std::vector<species>, int);

        void run_simulation(field *, field *, field *, field *, std::vector<species>, std::string);

        // MPI variables
        int grid_rank, rank;                              // rank of the current proces in the virtual grid
        int grid_top, grid_bottom, grid_left, grid_right; // ranks of the neighbouring processes
        int grid_ne, grid_se, grid_nw, grid_sw;           // ranks of diagonal processes: NE, SE, NW, SW

        int Nspecies;
        std::vector<int> Npart;
        std::vector<double> charge, mass, temp, vxfluid, vyfluid;

        // HDF5 variables
        hid_t file_field, dataset_field, dataspace_field;
        hid_t group_rank, dataset_rank, dataspace_rank;
        hid_t dataspace_part, dataset_part;
        hid_t group_charge, group_Ex, group_Ey;
        hid_t part_id;
        herr_t status_h5;
        hid_t group_creation_plist;
        std::vector<hid_t> h5_vec_group;

        int sim_true;

    private:
        MPI_Datatype exchange_field_type[2]; // MPI_datatype for exchange of buffer cell data
        MPI_Comm grid_comm;                  // grid COMMUNICATOR
        int offset[2];                       // offset for cell numbering for subdomains
        int wrap_around[2];
        MPI_Status status_mpi;

        // MPI_Datatype exchange_part_type;
        MPI_Aint offsets[7]; // it evaluates to the offset (in bytes) of a given member within a struct or union type
        const int nitems = 7;
        MPI_Datatype exchange_part_type;

        // Simulation variables
        double aspect, xlen;
        double *X_guard_data, *Y_guard_data;

        // Simulation time variables
        float time1, time2, total_time;
        float setup_time, hdf5_time;
        float particle_time, field_time; 
    };
}
#endif