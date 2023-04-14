//FCPIC - 2D Particle-in-Cell code using MPI
//Guilherme Crispim, João Palma, José Afonso
//Advanced Topics in Computational Physics, 2023, IST

//File simulation_mpi.cc:
//Implementation of all MPI communication related functions 
//in the class Simulation

#include "simulation.hh"
#include <fstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <unistd.h>

namespace FCPIC
{
    /* 
        Function setup_proc_grid
        + Sets a cartesian topology for the MPI processes, using
        existing MPI functions (MPI_Cart_...). This facilitates the
        organization
    */
    void simulation::setup_proc_grid()
    {

        int reorder = 1; //Reorder the rank IDs of the processes inside the grid

        //Creates the cartesian topology of the dimensions provided by the user
        MPI_Cart_create(MPI_COMM_WORLD, 2, grid, wrap_around, reorder, &grid_comm);

        //Retrieve new rank and cartesian coordinates of this process
        MPI_Comm_rank(grid_comm, &grid_rank);
        MPI_Cart_coords(grid_comm, grid_rank, 2, grid_coord);

        //Calculate ranks of neighboring processes in the grid
        MPI_Cart_shift(grid_comm, 0, 1, &grid_left, &grid_right);
        MPI_Cart_shift(grid_comm, 1, 1, &grid_bottom, &grid_top);

        //Datatype for MPI Communication in fields
        //Datatype for horizontal data exchange
        MPI_Type_vector(N_int_y, 1, N_x, MPI_FLOAT, &exchange_field_type[X_DIR]);
        MPI_Type_commit(&exchange_field_type[X_DIR]);

        //Datatype for vertical data exchange
        MPI_Type_vector(N_x, 1, 1, MPI_FLOAT, &exchange_field_type[Y_DIR]);
        MPI_Type_commit(&exchange_field_type[Y_DIR]);

        //Datatype for MPI particle communication
        int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
        MPI_Datatype types[7] = {MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_CXX_BOOL};

        offsets[0] = offsetof(part, ix);
        offsets[1] = offsetof(part, iy);
        offsets[2] = offsetof(part, x);
        offsets[3] = offsetof(part, y);
        offsets[4] = offsetof(part, ux);
        offsets[5] = offsetof(part, uy);
        offsets[6] = offsetof(part, flag);

        MPI_Type_create_struct(nitems, blocklengths, offsets, types, &exchange_part_type);
        MPI_Type_commit(&exchange_part_type);

        //Coordinates of the diagonal processes
        int coords_ne[2] = {grid_coord[X_DIR] + 1, grid_coord[Y_DIR] + 1};
        int coords_se[2] = {grid_coord[X_DIR] + 1, grid_coord[Y_DIR] - 1};
        int coords_nw[2] = {grid_coord[X_DIR] - 1, grid_coord[Y_DIR] + 1};
        int coords_sw[2] = {grid_coord[X_DIR] - 1, grid_coord[Y_DIR] - 1};

        get_diagonal_rank(coords_ne, grid_ne);
        get_diagonal_rank(coords_se, grid_se);
        get_diagonal_rank(coords_nw, grid_nw);
        get_diagonal_rank(coords_sw, grid_sw);
    }

    /* 
        Function get_diagonal_rank
        Inputs: coordinates of the diagonal process, id to be written 
        with the rank
        + Gets the rank of the diagonal processes to the grid
    */
    void simulation::get_diagonal_rank(int *coords, int &id_proc)
    {
        // both wrap around in X_DIR and Y_DIR must be 0 or 1 at the same time
        if (!wrap_around[X_DIR] && !wrap_around[Y_DIR] && (coords[0] >= grid[0] || coords[1] >= grid[1] || coords[0] < 0 || coords[1] < 0))
            id_proc = MPI_PROC_NULL;
        else
            MPI_Cart_rank(grid_comm, coords, &id_proc);
    }

    /* 
        Function exchange_phi_buffers
        Inputs: electric potential "field" object  
        + Exchanges the guard cells between adjacent processes 
        for the Jacobi iteration
    */
    void simulation::exchange_phi_buffers(field *phi)
    {
        //Indexes for fixing the first memory position of the guard rows
        int i, j;

        i = 1;
        j = 0;

        //Communication stream leftward
        MPI_Sendrecv(&phi->val[WEST_BOUND], 1, exchange_field_type[X_DIR], grid_left, 0,
                     &phi->val[EAST_GUARD], 1, exchange_field_type[X_DIR], grid_right, 0,
                     grid_comm, &status_mpi);

        //Communication stream rightward
        MPI_Sendrecv(&phi->val[EAST_BOUND], 1, exchange_field_type[X_DIR], grid_right, 0,
                     &phi->val[WEST_GUARD], 1, exchange_field_type[X_DIR], grid_left, 0,
                     grid_comm, &status_mpi);

        //Communication stream upward
        MPI_Sendrecv(&phi->val[NORTH_BOUND], 1, exchange_field_type[Y_DIR], grid_top, 0,
                     &phi->val[SOUTH_GUARD], 1, exchange_field_type[Y_DIR], grid_bottom, 0,
                     grid_comm, &status_mpi);

        //Communication stream downward
        MPI_Sendrecv(&phi->val[SOUTH_BOUND], 1, exchange_field_type[Y_DIR], grid_bottom, 0,
                     &phi->val[NORTH_GUARD], 1, exchange_field_type[Y_DIR], grid_top, 0,
                     grid_comm, &status_mpi);
    }

    /* 
        Function exchange_charge_buffers
        Inputs: charge "field" object  
        + Adds charge value deposited in guard cells
        of adjacent processes 
    */
    void simulation::exchange_charge_buffers(field *charge)
    {
        //Indexes for fixing the first memory position of the guard rows
        int i = 1;
        int j = 0;

        //Communication stream upward
        MPI_Sendrecv(&charge->val[NORTH_GUARD], 1, exchange_field_type[Y_DIR], grid_top, 0,
                     &Y_guard_data[0], N_x, MPI_FLOAT, grid_bottom, 0, grid_comm, &status_mpi);
        if (grid_bottom != MPI_PROC_NULL)
            charge->reduceSouthBound(Y_guard_data);

        //Communication stream downward
        MPI_Sendrecv(&charge->val[SOUTH_GUARD], 1, exchange_field_type[Y_DIR], grid_bottom, 0,
                     &Y_guard_data[0], N_x, MPI_FLOAT, grid_top, 0, grid_comm, &status_mpi);
        if (grid_top != MPI_PROC_NULL)
            charge->reduceNorthBound(Y_guard_data);

        //Communication stream leftward
        MPI_Sendrecv(&charge->val[WEST_GUARD], 1, exchange_field_type[X_DIR], grid_left, 0,
                     &X_guard_data[0], N_int_y, MPI_FLOAT, grid_right, 0, grid_comm, &status_mpi);
        if (grid_right != MPI_PROC_NULL)
            charge->reduceEastBound(X_guard_data);

        //Communication stream rightward
        MPI_Sendrecv(&charge->val[EAST_GUARD], 1, exchange_field_type[X_DIR], grid_right, 0,
                     &X_guard_data[0], N_int_y, MPI_FLOAT, grid_left, 0, grid_comm, &status_mpi);
        if (grid_left != MPI_PROC_NULL)
            charge->reduceWestBound(X_guard_data);
    }
    
    /* 
        Function exchange_particle_buffers
        Inputs: species object  
        + Sends out of bounds particles to the respective process
    */
    void simulation::exchange_particles_buffers(species *lepton)
    {
        lepton->size_recv_north = 0;
        lepton->size_recv_south = 0;
        lepton->size_recv_east = 0;
        lepton->size_recv_west = 0;

        lepton->size_recv_ne = 0;
        lepton->size_recv_se = 0;
        lepton->size_recv_nw = 0;
        lepton->size_recv_sw = 0;

        //Communication to determine the size of the arrays of each buffer
        MPI_Sendrecv(&(lepton->size_send_north), 1, MPI_INT, grid_top, 0, &(lepton->size_recv_south), 1, MPI_INT, grid_bottom, 0, grid_comm, &status_mpi);
        MPI_Sendrecv(&(lepton->size_send_south), 1, MPI_INT, grid_bottom, 0, &(lepton->size_recv_north), 1, MPI_INT, grid_top, 0, grid_comm, &status_mpi);
        MPI_Sendrecv(&(lepton->size_send_west), 1, MPI_INT, grid_left, 0, &(lepton->size_recv_east), 1, MPI_INT, grid_right, 0, grid_comm, &status_mpi);
        MPI_Sendrecv(&(lepton->size_send_east), 1, MPI_INT, grid_right, 0, &(lepton->size_recv_west), 1, MPI_INT, grid_left, 0, grid_comm, &status_mpi);

        MPI_Sendrecv(&(lepton->size_send_ne), 1, MPI_INT, grid_ne, 0, &(lepton->size_recv_sw), 1, MPI_INT, grid_sw, 0, grid_comm, &status_mpi);
        MPI_Sendrecv(&(lepton->size_send_sw), 1, MPI_INT, grid_sw, 0, &(lepton->size_recv_ne), 1, MPI_INT, grid_ne, 0, grid_comm, &status_mpi);
        MPI_Sendrecv(&(lepton->size_send_nw), 1, MPI_INT, grid_nw, 0, &(lepton->size_recv_se), 1, MPI_INT, grid_se, 0, grid_comm, &status_mpi);
        MPI_Sendrecv(&(lepton->size_send_se), 1, MPI_INT, grid_se, 0, &(lepton->size_recv_nw), 1, MPI_INT, grid_nw, 0, grid_comm, &status_mpi);

        //Allocate memory for the vectors that are going to receive the MPI particles
        part recv_dummy;
        recv_dummy.ix = -1; //Set to -1 as a way to check later if there was "actual" communication
        recv_dummy.iy = -1;
        lepton->recv_buffer_east.assign(lepton->size_recv_east, recv_dummy);
        lepton->recv_buffer_west.assign(lepton->size_recv_west, recv_dummy);
        lepton->recv_buffer_north.assign(lepton->size_recv_north, recv_dummy);
        lepton->recv_buffer_south.assign(lepton->size_recv_south, recv_dummy);

        lepton->recv_buffer_ne.assign(lepton->size_recv_ne, recv_dummy);
        lepton->recv_buffer_nw.assign(lepton->size_recv_nw, recv_dummy);
        lepton->recv_buffer_se.assign(lepton->size_recv_se, recv_dummy);
        lepton->recv_buffer_sw.assign(lepton->size_recv_sw, recv_dummy);

        //Buffers Communication
        //All traffic in direction "top"
        MPI_Sendrecv(&(lepton->send_buffer_north[0]), lepton->send_buffer_north.size(), exchange_part_type, grid_top, 0, &(lepton->recv_buffer_south[0]), lepton->size_recv_south, exchange_part_type, grid_bottom, 0, grid_comm, MPI_STATUS_IGNORE);
        //All traffic in direction "bottom"
        MPI_Sendrecv(&(lepton->send_buffer_south[0]), lepton->send_buffer_south.size(), exchange_part_type, grid_bottom, 0, &(lepton->recv_buffer_north[0]), lepton->size_recv_north, exchange_part_type, grid_top, 0, grid_comm, MPI_STATUS_IGNORE);
        //All traffic in direction "right"
        MPI_Sendrecv(&(lepton->send_buffer_west[0]), lepton->send_buffer_west.size(), exchange_part_type, grid_left, 0, &(lepton->recv_buffer_east[0]), lepton->size_recv_east, exchange_part_type, grid_right, 0, grid_comm, MPI_STATUS_IGNORE);
        //All traffic in direction "left"
        MPI_Sendrecv(&(lepton->send_buffer_east[0]), lepton->send_buffer_east.size(), exchange_part_type, grid_right, 0, &(lepton->recv_buffer_west[0]), lepton->size_recv_west, exchange_part_type, grid_left, 0, grid_comm, MPI_STATUS_IGNORE);

        //All traffic in direction "ne-sw"
        MPI_Sendrecv(&(lepton->send_buffer_ne[0]), lepton->send_buffer_ne.size(), exchange_part_type, grid_ne, 0, &(lepton->recv_buffer_sw[0]), lepton->size_recv_sw, exchange_part_type, grid_sw, 0, grid_comm, &status_mpi);
        //All traf\fic in direction "sw-ne"
        MPI_Sendrecv(&(lepton->send_buffer_sw[0]), lepton->send_buffer_sw.size(), exchange_part_type, grid_sw, 0, &(lepton->recv_buffer_ne[0]), lepton->size_recv_ne, exchange_part_type, grid_ne, 0, grid_comm, &status_mpi);
        //All traf\fic in direction "se-nw"
        MPI_Sendrecv(&(lepton->send_buffer_se[0]), lepton->send_buffer_se.size(), exchange_part_type, grid_se, 0, &(lepton->recv_buffer_nw[0]), lepton->size_recv_nw, exchange_part_type, grid_nw, 0, grid_comm, &status_mpi);
        //All traffic in direction "nw-se"
        MPI_Sendrecv(&(lepton->send_buffer_nw[0]), lepton->send_buffer_nw.size(), exchange_part_type, grid_nw, 0, &(lepton->recv_buffer_se[0]), lepton->size_recv_se, exchange_part_type, grid_se, 0, grid_comm, &status_mpi);
    }
}