#ifndef __FCPICbase__
#define __FCPICbase__

namespace FCPIC{
    class FCPIC_base{
        protected:

        FCPIC_base(){};

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

        float Lref, Vref, Tref, Nref, Tempref;
        float dx, dy;
        float dt, simtime;
        int N_total_x, N_total_y;
        int N_int_x, N_int_y;
        int N_x, N_y, N;
        
        int bc;

        //MPI
        int n_Procs;
        int grid_coord[2];
        int grid[2];
    };
}

#endif