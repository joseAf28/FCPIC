#ifndef _SIMULATION_
#define _SIMULATION_

//FCPIC - 2D Particle-in-Cell code using MPI
//Guilherme Crispim, João Palma, José Afonso
//Advanced Topics in Computational Physics, 2023, IST

//File simulation.hh:
//Declaration of class simulation

#include "mpi.h"
#include "hdf5.h"
#include "species.hh"
#include "field.hh"
#include "FCPIC_base.hh"

namespace FCPIC
{
    //Class simulation
    //Includes all numerical methods for evolving the simulation variables,
    //IO interaction with the user, all input handling and setting of the 
    //parameters, all MPI communication and writing of the data to HDF5 files,
    //for diagnostics. It derives from an FCPIC_base object, which is empty and
    //written with all parameters of the simulation. This FCPIC_base object is
    //thus sent to all other simulation objects, for consistency of all elements
    class simulation : public FCPIC_base
    {
    public:
        // Constructor/destructor
        simulation(int, char **);
        ~simulation() override;

        //IO handling
        void readArgs(int, char **);
        void getParamsfromFile(std::string, std::vector<bool> *);
        void setParams();
        void confirmParams();
        void printHelp();
        std::string print_SI(float, int);
        void printTitle();
        void printProgress(float);
        void printTime(std::string);

        //Time functions
        void setTime(float &);
        void setTime();

        //MPI functions
        void setup_proc_grid();
        void get_diagonal_rank(int *, int &);

        void exchange_phi_buffers(field *);
        void exchange_charge_buffers(field *);
        void exchange_particles_buffers(species *);

        //Numerical methods
        void get_charge(field *, species *);
        void jacobi(field *, field *);

        void set_E_value(field *, field *, field *);

        void field_interpolate(field *, field *, float &, float &, part *);
        void init_pusher(field *, field *, species *);
        void particle_pusher(field *, field *, species *);

        //HDF5 functions
        void setupHDF5(std::string);
        void closeHDF5(std::string);
        void writeChargeHDF5(field *, int);
        void writeExHDF5(field *, int);
        void writeEyHDF5(field *, int);
        void writePartHDF5(std::vector<species> &, int);

        //Function progressing the simulation
        void run_simulation(field *, field *, field *, field *, std::vector<species> &, std::string);

        //MPI variables
        int grid_rank, rank;                              // rank of the current proces in the virtual grid
        int grid_top, grid_bottom, grid_left, grid_right; // ranks of the neighbouring processes
        int grid_ne, grid_se, grid_nw, grid_sw;           // ranks of diagonal processes: NE, SE, NW, SW

        int Nspecies;
        std::vector<int> Npart, Nypart, Nxpart, rand_true;
        std::vector<float> charge, mass, temp, vxfluid, vyfluid;

        //HDF5 variables
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
        float aspect, xlen;
        float *X_guard_data, *Y_guard_data;

        // Simulation time variables
        float time1, time2, total_time;
        float setup_time, hdf5_time;
        float particle_time, field_time;
    };
}
#endif