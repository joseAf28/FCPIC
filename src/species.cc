//FCPIC - 2D Particle-in-Cell code using MPI
//Guilherme Crispim, João Palma, José Afonso
//Advanced Topics in Computational Physics, 2023, IST

//File field.cc:
//Implementation of all methods of class Species

#include <random>
#include <algorithm>
#include <iomanip>
#include "species.hh"
#include "math.h"
#include "mpi.h"

namespace FCPIC
{   

    /* 
        Constructor of class species
        Inputs: charge, mass, temperature, array for components of 
        fluid velocity (all normalized), array for number of particles 
        per X and Y direction, and FCPIC_base object with all parameters 
        previously set
        + Creates an initial configuration for a species disposed
        in a uniform grid of ppd[0] (along x) x ppd[1] (along y)
        in the simulation box, and with normal distribution of
        velocities with avg = vf and sigma = v_thermal. Each call of this
        constructor (once per process) only handles the particles in the 
        current process, avoiding comms
    */
    species::species(float charge, float mass, float temp, float *vf, int *ppd, FCPIC_base const *base) : FCPIC_base(*base), q(charge), m(mass)
    {
        //Number of particles in the entire simulation
        np_sim = ppd[0] * ppd[1];
        
        //Distance between consecutive particles in both directions
        float d_part_x = (dx * ((float)N_total_x)) / ((float)ppd[0]);
        float d_part_y = (dy * ((float)N_total_y)) / ((float)ppd[1]);

        //Aux vectors for pushing back X and Y positions
        std::vector<float> x_aux, y_aux;
        x_aux.clear();
        y_aux.clear();

        //Distinction between BCs. Done due to the extra physical row of
        //physical cells in conductive, requires an extra offset of dx
        if (bc == PERIODIC)
        {
            float local_val;
            float L_int_x = ((float)N_int_x) * dx; //Physical lengths of
            float L_int_y = ((float)N_int_y) * dy; //a process domain

            //Testing if particle grid positions along X land in the domain
            //of the current process
            for (int i = 0; i < ppd[0]; i++)
            {
                local_val = (((float)i) + .5) * d_part_x - ((float)grid_coord[0]) * L_int_x;
                if (local_val >= 0. && local_val < L_int_x)
                    x_aux.push_back(local_val);
            }

            //Testing if particle grid positions along Y land in the domain
            //of the current process
            for (int i = 0; i < ppd[1]; i++)
            {
                local_val = (((float)i) + .5) * d_part_y - ((float)grid_coord[1]) * L_int_y;
                if (local_val >= 0. && local_val < L_int_y)
                    y_aux.push_back(local_val);
            }
        }
        //Conductive BC case
        else
        {
            float local_val;
            float L_int_x = ((float)N_int_x) * dx;
            float L_int_y = ((float)N_int_y) * dy;

            for (int i = 0; i < ppd[0]; i++)
            {
                //Subtracting an extra dx
                local_val = (((float)i) + .5) * d_part_x - ((float)grid_coord[0]) * L_int_x - dx;
                if (local_val >= 0. && local_val < L_int_x)
                    x_aux.push_back(local_val);
                else if (grid_coord[0] == 0 && local_val < 0.)
                    x_aux.push_back(local_val);
            }

            for (int i = 0; i < ppd[1]; i++)
            {
                //Subtracting an extra dy
                local_val = (((float)i) + .5) * d_part_y - ((float)grid_coord[1]) * L_int_y - dy;
                if (local_val >= 0. && local_val < L_int_y)
                    y_aux.push_back(local_val);
                else if (grid_coord[1] == 0 && local_val < 0.)
                    y_aux.push_back(local_val);
            }
        }

        //Initial number of particles in the current process domain
        np = x_aux.size() * y_aux.size();

        //Reserve space for the arrays of particles
        part A;
        A.flag = BULK;
        vec.reserve(3 * np); //Assumption for the minimum reserved space
        vec.assign(np, A);

        //Aux vectors for communication of particles
        send_buffer_north.reserve(np); 
        send_buffer_south.reserve(np); 
        send_buffer_east.reserve(np);
        send_buffer_west.reserve(np);

        send_buffer_ne.reserve(np);
        send_buffer_se.reserve(np); 
        send_buffer_nw.reserve(np);
        send_buffer_sw.reserve(np);

        recv_buffer_north.reserve(np);
        recv_buffer_south.reserve(np);
        recv_buffer_east.reserve(np);
        recv_buffer_west.reserve(np);

        recv_buffer_ne.reserve(np);
        recv_buffer_se.reserve(np);
        recv_buffer_nw.reserve(np);
        recv_buffer_sw.reserve(np);

        //Random number generator for velocity generation
        std::random_device dev;
        std::mt19937_64 rng(dev());
        std::normal_distribution<float> rand_gauss(0, 1);

        //Thermal velocity normalized
        float vth = sqrt(temp/mass);

        int n_a = 0; //Aux counter of particles

        //Ranging over the selected domain
        for (int i = 0; i < y_aux.size(); i++)
            for (int j = 0; j < x_aux.size(); j++)
            {
                //The position is stored in x and y, in reference to the
                //point of index (1,1) in the current process. Given x and y
                //are not necessarily in the interval [0,dx[ and [0,dy[,
                //respectively, the positions need to be handled by the 
                //function advance_cell()
                vec.at(n_a).x = x_aux[j];
                vec.at(n_a).y = y_aux[i];
                vec.at(n_a).ix = 1;
                vec.at(n_a).iy = 1;
                vec.at(n_a).ux = vf[0] + vth * rand_gauss(rng);
                vec.at(n_a).uy = vf[1] + vth * rand_gauss(rng);
                n_a += 1;
            }

        x_aux.clear();
        y_aux.clear();
    }

    /* 
        Constructor of class species
        Inputs: charge, mass, temperature, array for components of 
        fluid velocity (all normalized), array for number of particles 
        per X and Y direction, and FCPIC_base object with all parameters 
        previously set
        + Creates an initial configuration for a species disposed
        with random positions with uniform distribution, and with 
        normal distribution of velocities with avg = vf and 
        sigma = v_thermal. Each call of this object (once per process),
        will generate [total number of particles]/[number of processes]
        particles at any position of the entire simulation. Hence,
        communications are necessary right after this.
    */
    species::species(float charge, float mass, float temp, float *vf, int n_part, FCPIC_base const *base) : FCPIC_base(*base), q(charge), m(mass), np_sim(n_part)
    {
        
        //Initial number of particles in this current process
        np = np_sim / n_Procs;

        //Reserve space for the arrays of particles
        part A;
        A.flag = BULK;
        vec.reserve(np);
        vec.assign(np, A);

        //Aux vectors for communication of particles
        send_buffer_north.reserve(np); 
        send_buffer_south.reserve(np); 
        send_buffer_east.reserve(np);
        send_buffer_west.reserve(np);

        send_buffer_ne.reserve(np);
        send_buffer_se.reserve(np); 
        send_buffer_nw.reserve(np);
        send_buffer_sw.reserve(np);

        recv_buffer_north.reserve(np);
        recv_buffer_south.reserve(np);
        recv_buffer_east.reserve(np);
        recv_buffer_west.reserve(np);

        recv_buffer_ne.reserve(np);
        recv_buffer_se.reserve(np);
        recv_buffer_nw.reserve(np);
        recv_buffer_sw.reserve(np);

        //Random number generator for position and velocity generation
        std::random_device dev;
        std::mt19937_64 rng(dev());
        std::normal_distribution<float> rand_gauss(0, 1);
        std::uniform_real_distribution<float> rand_uniform(0.0, 1.0);

        //Thermal velocity normalized
        float vth = sqrt(temp/mass);

        
        if (bc == CONDUCTIVE)
        {
            for (int i = 0; i < np; i++)
            {
                //The position is stored in x and y, in reference to the
                //point of index (1,1) in the current process. Given x and y
                //are not necessarily in the interval [0,dx[ and [0,dy[,
                //respectively, the positions need to be handled by the 
                //function advance_cell()
                vec.at(i).x = (N_total_x * rand_uniform(rng) - ((float)grid_coord[0] * (float)N_int_x + 1.)) * dx; //Subtracting an extra dx for conductive BCs
                vec.at(i).y = (N_total_y * rand_uniform(rng) - ((float)grid_coord[1] * (float)N_int_y + 1.)) * dy; //Subtracting an extra dy for conductive BCs
                vec.at(i).ix = 1;
                vec.at(i).iy = 1;
                vec.at(i).ux = vf[0] + vth * rand_gauss(rng);
                vec.at(i).uy = vf[1] + vth * rand_gauss(rng);
            }
        }
        else
        {
            for (int i = 0; i < np; i++)
            {
                vec.at(i).x = (N_total_x * rand_uniform(rng) - ((float)grid_coord[0] * (float)N_int_x)) * dx;
                vec.at(i).y = (N_total_y * rand_uniform(rng) - ((float)grid_coord[1] * (float)N_int_y)) * dy;
                vec.at(i).ix = 1;
                vec.at(i).iy = 1;
                vec.at(i).ux = vf[0] + vth * rand_gauss(rng);
                vec.at(i).uy = vf[1] + vth * rand_gauss(rng);
            }
        }
    }

    /* 
        Destructor of class species
        + Free all memory for the particles stored and for
        the aux buffers
    */
    species::~species()
    {
        vec.clear();

        send_buffer_north.clear();
        send_buffer_south.clear();
        send_buffer_east.clear();
        send_buffer_west.clear();

        send_buffer_nw.clear();
        send_buffer_sw.clear();
        send_buffer_ne.clear();
        send_buffer_se.clear();

        recv_buffer_north.clear();
        recv_buffer_south.clear();
        recv_buffer_east.clear();
        recv_buffer_west.clear();

        recv_buffer_ne.clear();
        recv_buffer_se.clear();
        recv_buffer_sw.clear();
        recv_buffer_nw.clear();
    }

    /* 
        Function advance_cell
        Inputs: ranks of current and adjacent processes
        ranks_mpi[0] - curr. rank, ranks_mpi[1] - top, 
        ranks_mpi[2] - bottom, ranks_mpi[3] - right, 
        ranks_mpi[4] - left
        Output: count of the number of processes with particles
        ready for being sent. When 0, we are sure we don't need to
        call this function again
        + Given initial generation of particles and 
        particle advancing methods only change the x and
        y entries of each particle, this function reduces
        both variables to the [0,dx[ and [0,dy[ interval,
        with index of nearest bottom left grid point ix and
        iy.
        + After reducing, this function checks and operates
        reflections of particles when reaching physical boundaries.
        + It also checks if particles need to be sent to other processes,
        and heads them to the respective send buffer. This operation is only 
        followed between adjacent processes (including diagonals). However,
        this function can be iterated the number of times necessary, until
        the message passes from process to process to the desired point. 
        This function, followed by the respective communication functions,
        must be iterated until the output flag is 0.
    */
    int species::advance_cell(int *ranks_mpi)
    { 
        int changes_made = 0; //Local flag testing if communication is necessary
        int global_changes_made = 0; //Global flag with the number of procs requiring comms
        //Cycle over all local particles
        for (int counter = 0; counter < np; counter++)
        {
            //If these conditions are met, the particle is already reduced
            //and in the correct process, so it can be skipped
            if (vec[counter].x >= 0 && vec[counter].x < dx &&
                vec[counter].y >= 0 && vec[counter].y < dy &&
                vec[counter].ix >= 0 && vec[counter].ix <= N_int_x &&
                vec[counter].iy >= 0 && vec[counter].iy <= N_int_y)
                continue;

            //Reduction of x variable
            vec[counter].ix += floor(vec[counter].x / dx);
            vec[counter].x = fmod(vec[counter].x, dx);

            if (vec[counter].x < 0)
                vec[counter].x += dx;

            //Reduction of y variable
            vec[counter].iy += floor(vec[counter].y / dy);
            vec[counter].y = fmod(vec[counter].y, dy);

            if (vec[counter].y < 0)
                vec[counter].y += dy;

            //Testing if the process has a physical boundary
            //in the respective direction (sending vs reflecting)
            bool flag_top = ranks_mpi[1] == MPI_PROC_NULL;
            bool flag_bottom = ranks_mpi[2] == MPI_PROC_NULL;
            bool flag_right = ranks_mpi[3] == MPI_PROC_NULL;
            bool flag_left = ranks_mpi[4] == MPI_PROC_NULL;

            //Testing if particle is out of bounds
            bool ixmin_cond = vec[counter].ix <= -1;
            bool ixmax_cond = vec[counter].ix > N_int_x;
            bool iymin_cond = vec[counter].iy <= -1;
            bool iymax_cond = vec[counter].iy > N_int_y;

            bool send_N_true = false;
            bool send_S_true = false;
            bool send_W_true = false;
            bool send_E_true = false;

            if (ixmin_cond) //Left 
            {
                if (flag_left) //Reflection
                {
                    vec[counter].x = dx - vec[counter].x;
                    vec[counter].ix = -vec[counter].ix - 1;
                    vec[counter].ux = -vec[counter].ux;
                }
                else //Communication
                {
                    vec[counter].ix = vec[counter].ix + N_int_x;
                    send_W_true = true;
                    changes_made = true;
                }
            }

            if (iymin_cond) //Bottom
            {
                if (flag_bottom) //Reflection
                {
                    vec[counter].y = dy - vec[counter].y;
                    vec[counter].iy = -vec[counter].iy - 1;
                    vec[counter].uy = -vec[counter].uy;
                }
                else //Communication
                {
                    vec[counter].iy = vec[counter].iy + N_int_y;
                    send_S_true = true;
                    changes_made = 1;
                }
            }

            if (ixmax_cond) //Right
            {
                if (flag_right) //Reflection
                {
                    vec[counter].x = dx - vec[counter].x;
                    vec[counter].ix = 2 * N_int_x - vec[counter].ix + 1;
                    vec[counter].ux = -vec[counter].ux;
                }
                else //Communication
                {
                    vec[counter].ix = vec[counter].ix - N_int_x;
                    send_E_true = true;
                    changes_made = 1;
                }
            }

            if (iymax_cond) //Top
            {
                if (flag_top) //Reflection
                {
                    vec[counter].y = dy - vec[counter].y;
                    vec[counter].iy = 2 * N_int_y - vec[counter].iy + 1;
                    vec[counter].uy = -vec[counter].uy;
                }
                else //Communication
                {
                    vec[counter].iy = vec[counter].iy - N_int_y;
                    send_N_true = true;
                    changes_made = 1;
                }
            }

            //Adding particles to respective sending buffer
            //N buffer
            if ((!send_W_true) && (!send_E_true) && send_N_true)
            {
                send_buffer_north.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            //NE buffer
            else if (send_N_true && send_E_true)
            {
                send_buffer_ne.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            //NW buffer
            else if (send_N_true && send_W_true)
            {
                send_buffer_nw.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            //S buffer
            else if ((!send_W_true) && (!send_E_true) && send_S_true)
            {
                send_buffer_south.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            //SE buffer
            else if (send_S_true && send_E_true)
            {
                send_buffer_se.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            //SW buffer
            else if (send_S_true && send_W_true)
            {
                send_buffer_sw.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            //E buffer
            else if (send_E_true && (!send_N_true) && (!send_S_true))
            {
                send_buffer_east.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            //W buffer
            else if (send_W_true && (!send_N_true) && (!send_S_true))
            {
                send_buffer_west.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            else{}
        }
        //Communication to update global flag. Even if this process does not need to send anything, 
        //it may have to receive particles. Without this flag, it would not predict and communication
        //and a deadlock would be met
        MPI_Allreduce(&changes_made, &global_changes_made, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        return global_changes_made;
    }

    /* 
        Function prepare_buffer
        + Erases all particles being communicated from the local list
    */
    void species::prepare_buffer()
    {
        //Determine the size of the arrays that Exchange particles in MPI
        size_send_north = send_buffer_north.size();
        size_send_south = send_buffer_south.size();
        size_send_east = send_buffer_east.size();
        size_send_west = send_buffer_west.size();

        size_send_ne = send_buffer_ne.size();
        size_send_se = send_buffer_se.size();
        size_send_nw = send_buffer_nw.size();
        size_send_sw = send_buffer_sw.size();

        //Delete all particles that were sent to exchange array buffers
        vec.erase(std::remove_if(vec.begin(), vec.end(), [this](const part obj)
                                 { return (obj.flag == SEND); }),
                  vec.end());
    }

    /* 
        Function update_part_list
        + Adds all particles received from MPI comms to the local list and
        clears the aux buffers
    */
    void species::update_part_list()
    {
        //Include the new particles after the MPI data exchange
        for (int i = 0; i < recv_buffer_north.size(); i++)
        {
            vec.push_back(recv_buffer_north[i]);
        }
        for (int i = 0; i < recv_buffer_south.size(); i++)
        {
            vec.push_back(recv_buffer_south[i]);
        }
        for (int i = 0; i < recv_buffer_east.size(); i++)
        {
            vec.push_back(recv_buffer_east[i]);
        }
        for (int i = 0; i < recv_buffer_west.size(); i++)
        {
            vec.push_back(recv_buffer_west[i]);
        }
        for (int i = 0; i < recv_buffer_ne.size(); i++)
        {
            vec.push_back(recv_buffer_ne[i]);
        }
        for (int i = 0; i < recv_buffer_nw.size(); i++)
        {
            vec.push_back(recv_buffer_nw[i]);
        }
        for (int i = 0; i < recv_buffer_sw.size(); i++)
        {
            vec.push_back(recv_buffer_sw[i]);
        }
        for (int i = 0; i < recv_buffer_se.size(); i++)
        {
            vec.push_back(recv_buffer_se[i]);
        }

        //Update the particles number in the simulation after MPI exchange
        np = vec.size();

        //Clean old data from buffers
        recv_buffer_north.clear();
        recv_buffer_south.clear();
        recv_buffer_east.clear();
        recv_buffer_west.clear();

        recv_buffer_ne.clear();
        recv_buffer_se.clear();
        recv_buffer_nw.clear();
        recv_buffer_sw.clear();

        send_buffer_north.clear();
        send_buffer_south.clear();
        send_buffer_east.clear();
        send_buffer_west.clear();

        send_buffer_ne.clear();
        send_buffer_se.clear();
        send_buffer_nw.clear();
        send_buffer_sw.clear();
    }
}