#ifndef __FCPICbase__
#define __FCPICbase__

//FCPIC - 2D Particle-in-Cell code using MPI
//Guilherme Crispim, João Palma, José Afonso
//Advanced Topics in Computational Physics, 2023, IST

//File FCPIC_base.hh:
//Declaration of class FCPIC_base

namespace FCPIC{

    //Class FCPIC_base
    //Includes all important simulation parameters that need to be accessed by all
    //parts of the simulation. All classes of the project need to derive from this
    class FCPIC_base{
        protected:

        //Default constructor
        //To be called in the derived class that is going to fill in the parameters
        //in the first place
        FCPIC_base(){};

        //Copy constructor
        //To be called by derived classes once the parameters are filled in
        FCPIC_base(FCPIC_base const & obj) : Lref(obj.Lref),
        Vref(obj.Vref), Tref(obj.Tref), Nref(obj.Nref), Tempref(obj.Tempref),
        N_total_x(obj.N_total_x), N_total_y(obj.N_total_y),
        N_int_x(obj.N_int_x), N_int_y(obj.N_int_y),
        N_x(obj.N_x), N_y(obj.N_y), N(obj.N),
        dx(obj.dx), dy(obj.dy), dt(obj.dt), 
        n_Procs(obj.n_Procs), simtime(obj.simtime), bc(obj.bc){

            grid_coord[0] = obj.grid_coord[0];
            grid_coord[1] = obj.grid_coord[1];

            grid[0] = obj.grid[0];
            grid[1] = obj.grid[1];
        };

        virtual ~FCPIC_base(){};

        float Lref, Vref, Tref, Nref, Tempref; //Normalizations
        float dx, dy;                          //Spatial discretization
        float dt, simtime;                     //Time discretization and sim time
        int N_total_x, N_total_y;              //Number of total grid cells per direction
        int N_int_x, N_int_y;                  //Number of physical grid cells per dir. per proc.
        int N_x, N_y, N;                       //N_int with the 2 added guard cells per dir.
        int bc;                                //Boundary condition type

        int n_Procs;                           //Number of MPI procs
        int grid_coord[2];                     //Coordinates of the current process
        int grid[2];                           //Dimensions of the MPI process grid
    };
}

#endif