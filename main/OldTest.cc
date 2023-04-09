#include "species.hh"
#include "mpi.h"
#include "simulation.hh"

// Direction
typedef enum
{
    X_DIR,
    Y_DIR,
} Direction;

int main(int argc, char **argv)
{
    ///!MPI COMMUNICATION**************************
    MPI_Init(&argc, &argv);

    // MPI variables
    int proc_rank;                                            // rank of the current process
    int P_grid_rank;                                          // rank of the current proces in the virtual grid
    int P_grid_top, P_grid_bottom, P_grid_left, P_grid_right; // ranks of the neighbouring processes
    int n_Procs;                                              // total number of processes
    int P_grid[2];                                            // virtual grid dimensions
    MPI_Comm grid_comm;                                       // grid COMMUNICATOR
    int offset[2];                                            // offset for cell numbering for subdomains
    int P_grid_coord[2];                                      // coordinates of the process in the virtual grid
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    int wrap_around[2];
    int reorder;

    // retrieve the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &n_Procs);

    // number of processes per row and column
    P_grid[X_DIR] = 2;
    P_grid[Y_DIR] = 2;

    if (P_grid[X_DIR] * P_grid[Y_DIR] != n_Procs)
        std::cout << "Error: Mismatch of number of processes and process grid" << std::endl;

    // creating a virtual process grid topology
    wrap_around[X_DIR] = 0;
    wrap_around[Y_DIR] = 0; // dont connect first and last processes
    reorder = 1;            // reorder process ranks

    // creates a new communicator, grid_comm
    MPI_Cart_create(MPI_COMM_WORLD, 2, P_grid, wrap_around, reorder, &grid_comm);

    // retrieve new rank and cartesian coordinates of this process
    MPI_Comm_rank(grid_comm, &P_grid_rank);
    MPI_Cart_coords(grid_comm, P_grid_rank, 2, P_grid_coord);

    // calculate ranks of neighboring processes in the grid
    MPI_Cart_shift(grid_comm, 1, 1, &P_grid_left, &P_grid_right);
    MPI_Cart_shift(grid_comm, 0, 1, &P_grid_bottom, &P_grid_top);

    std::cout << "rank: " << P_grid_rank << " right: " << P_grid_right << " left: " << P_grid_left << " top: " << P_grid_top << " bottom: " << P_grid_bottom << std::endl;
    std::cout << "Proc_null: " << MPI_PROC_NULL << std::endl;

    // Create Structure Data type MPI
    const int nitems = 8;
    int blocklengths[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[8] = {MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_CXX_BOOL};
    MPI_Datatype mpi_part;
    MPI_Aint offsets[8];

    offsets[0] = offsetof(part, ix); // it evaqluatez to the offset (in bytes) of a given member within a struct or union type
    offsets[1] = offsetof(part, iy);
    offsets[2] = offsetof(part, x);
    offsets[3] = offsetof(part, y);
    offsets[4] = offsetof(part, ux);
    offsets[5] = offsetof(part, uy);
    offsets[6] = offsetof(part, uz);
    offsets[7] = offsetof(part, flag);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_part);
    MPI_Type_commit(&mpi_part);

    ///!MPI COMMUNICATION**************************

    // Initialize the Species Class
    std::string name = "electron";

    int ppc[2] = {1, 1};
    int range[2] = {5, 5}; // number of cells in each direction

    float *vf = new float[3];
    float vth[3] = {0.0, 0.0, 0.0};

    vf[0] = -1.;
    vf[1] = -1.;
    vf[2] = -1.;

    // differentiate vectors
    if (proc_rank == 0) // 0
    {
        vf[0] = -1.;
        vf[1] = 0;
        vf[2] = 0;
    }
    if (proc_rank == 1) // 1
    {
        vf[0] = -1;
        vf[1] = -1;
        vf[2] = 0.;
    }
    if (proc_rank == 2) // 2
    {
        vf[0] = -1;
        vf[1] = -1;
        vf[2] = 0.;
    }
    if (proc_rank == 3) // 3
    {
        vf[0] = -1;
        vf[1] = -1;
        vf[2] = 0.;
    }

    FCPIC::field *Ex = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *Ey = new FCPIC::field(range[0] + 1, range[1] + 1);

    FCPIC::field *charge = new FCPIC::field(range[0] + 1, range[1] + 1); // intialize to zero in all entries
    FCPIC::field *phi = new FCPIC::field(range[0] + 1, range[1] + 1);

    species test(name, ppc, range, vf, vth, 1);
    if (proc_rank == 0)
    {
        for (int i = 1; i < test.vec.size(); i++)
            test.vec[i].flag = SEND;

        test.vec.erase(std::remove_if(test.vec.begin(), test.vec.end(), [&test](const part obj)
                                      { return (obj.flag == SEND); }),
                       test.vec.end());

        std::cout << proc_rank << "vec.size: " << test.vec.size() << std::endl;
    }

    if (proc_rank != 0)
    {
        for (int i = 0; i < test.vec.size(); i++)
            test.vec[i].flag = SEND;

        test.vec.erase(std::remove_if(test.vec.begin(), test.vec.end(), [&test](const part obj)
                                      { return (obj.flag == SEND); }),
                       test.vec.end());
        std::cout << proc_rank << "vec.size: " << test.vec.size() << std::endl;
    }
    test.set_x();
    test.set_u();

    // test.write_output_vec(0, P_grid_rank);

    // sim->set_E_value(phi, Ex, Ey);
    // //
    for (int i = 0; i < 500; i++)
    {
        int flags_coords_mpi[5] = {P_grid_rank, P_grid_top, P_grid_bottom, P_grid_right, P_grid_left};

        test.update_part_list();
        test.write_output_vec(P_grid_rank, i, P_grid_rank);
        test.init_pusher(Ex, Ey);
        test.particle_pusher(Ex, Ey);
        test.advance_cell(flags_coords_mpi);
        std::cout << "+++*+***+**+**+*+******++*+*+++++++++" << std::endl;

        //
        test.prepare_buffer();

        //!!! Size of the arrays to send
        MPI_Sendrecv(&(test.size_send_north), 1, MPI_INT, P_grid_top, 0, &(test.size_recv_south), 1, MPI_INT, P_grid_bottom, 0, grid_comm, &status);
        MPI_Sendrecv(&(test.size_send_south), 1, MPI_INT, P_grid_bottom, 0, &(test.size_recv_north), 1, MPI_INT, P_grid_top, 0, grid_comm, &status);
        MPI_Sendrecv(&(test.size_send_west), 1, MPI_INT, P_grid_left, 0, &(test.size_recv_east), 1, MPI_INT, P_grid_right, 0, grid_comm, &status);
        MPI_Sendrecv(&(test.size_send_east), 1, MPI_INT, P_grid_right, 0, &(test.size_recv_west), 1, MPI_INT, P_grid_left, 0, grid_comm, &status);

        part recv_dummy;
        recv_dummy.ix = -1;
        recv_dummy.iy = -1;

        test.recv_buffer_east.assign(test.size_recv_east, recv_dummy);
        test.recv_buffer_west.assign(test.size_recv_west, recv_dummy);
        test.recv_buffer_north.assign(test.size_recv_north, recv_dummy);
        test.recv_buffer_south.assign(test.size_recv_south, recv_dummy);

        // test.write_input_buffer(i, P_grid_rank);

        //! Buffers Communication
        MPI_Sendrecv(&(test.send_buffer_north[0]), test.send_buffer_north.size(), mpi_part, P_grid_top, 0, &(test.recv_buffer_south[0]), test.size_recv_south, mpi_part, P_grid_bottom, 0, grid_comm, &status);

        // All traf\fic in direction "top"
        MPI_Sendrecv(&(test.send_buffer_south[0]), test.send_buffer_south.size(), mpi_part, P_grid_bottom, 0, &(test.recv_buffer_north[0]), test.size_recv_north, mpi_part, P_grid_top, 0, grid_comm, &status);

        MPI_Sendrecv(&(test.send_buffer_west[0]), test.send_buffer_west.size(), mpi_part, P_grid_left, 0, &(test.recv_buffer_east[0]), test.size_recv_east, mpi_part, P_grid_right, 0, grid_comm, &status);

        // All traffic in direction "right"
        MPI_Sendrecv(&(test.send_buffer_east[0]), test.send_buffer_east.size(), mpi_part, P_grid_right, 0, &(test.recv_buffer_west[0]), test.size_recv_west, mpi_part, P_grid_left, 0, grid_comm, &status);

        std::cout << "size north: " << test.send_buffer_north.size() << std::endl;

        // test.write_input_buffer(8, P_grid_rank);
        // test.update_part_list();

        charge->setValue(0.f);

        test.get_charge(charge);
    }

    MPI_Finalize();

    return 0;
}