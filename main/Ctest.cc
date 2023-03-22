#include "species.hh"
#include "mpi.h"

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
    const int nitems = 7;
    int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[7] = {MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_part;
    MPI_Aint offsets[7];

    offsets[0] = offsetof(part, ix); // it evaqluatez to the offset (in bytes) of a given member within a struct or union type
    offsets[1] = offsetof(part, iy);
    offsets[2] = offsetof(part, x);
    offsets[3] = offsetof(part, y);
    offsets[4] = offsetof(part, ux);
    offsets[5] = offsetof(part, uy);
    offsets[6] = offsetof(part, uz);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_part);
    MPI_Type_commit(&mpi_part);
    ///!MPI COMMUNICATION**************************

    // Initialize the Soecies Class
    std::string name = "electron";

    int ppc[2] = {2, 2};
    int range[2] = {10, 10}; // number of cells in each direction

    float *vf = new float[3];
    float vth[3] = {0.0, 0.0, 0.0};

    // differentiate vectors
    if (proc_rank == 0) // 0
    {
        vf[0] = 0;
        vf[1] = 0;
        vf[2] = 0;
    }
    if (proc_rank == 1) // 1
    {
        vf[0] = 1;
        vf[1] = 1;
        vf[2] = 1;
    }
    if (proc_rank == 2) // 2
    {
        vf[0] = 2;
        vf[1] = 2;
        vf[2] = 2;
    }
    if (proc_rank == 3) // 3
    {
        vf[0] = 3;
        vf[1] = 3;
        vf[2] = 3;
    }

    species test(name, ppc, range, vf, vth);
    test.set_x();
    test.set_u();
    test.get_charge();
    test.write_output_vec(2, P_grid_rank);
    test.to_buffer();

    if (P_grid_rank == 2)
    {
        std::cout << "north buffer" << std::endl;
        for (int i = 0; i < test.send_buffer_north.size(); i++)
            std::cout << "i: " << i << " ix: " << test.send_buffer_north[i].ix << " iy: " << test.send_buffer_north[i].iy << std::endl;
    }
    // int len_mpi = test.send_buffer_north.size();
    int west_len_mpi = test.buffer_west_len;
    int east_len_mpi = test.buffer_east_len;
    int north_len_mpi = test.buffer_north_len;
    int south_len_mpi = test.buffer_south_len;

    int vert_len_mpi = std::max(north_len_mpi, south_len_mpi);
    int hor_len_mpi = std::max(east_len_mpi, west_len_mpi);

    std::cout << "north: " << north_len_mpi << " south: " << south_len_mpi << " east: " << east_len_mpi << " west: " << west_len_mpi << std::endl;

    //!! Send Data among ranks
    // All traffic in direction "top"
    MPI_Sendrecv(&(test.send_buffer_north[0]), north_len_mpi, mpi_part, P_grid_top, 0, &(test.recv_buffer_south[0]), vert_len_mpi, mpi_part, P_grid_bottom, 0, grid_comm, &status);

    // All traffic in direction "top"
    MPI_Sendrecv(&(test.send_buffer_south[0]), south_len_mpi, mpi_part, P_grid_bottom, 0, &(test.recv_buffer_north[0]), vert_len_mpi, mpi_part, P_grid_top, 0, grid_comm, &status);

    // All traffic in direction "left"
    MPI_Sendrecv(&(test.send_buffer_west[0]), west_len_mpi, mpi_part, P_grid_left, 0, &(test.recv_buffer_east[0]), hor_len_mpi, mpi_part, P_grid_right, 0, grid_comm, &status);

    // All traffic in direction "right"
    MPI_Sendrecv(&(test.send_buffer_east[0]), east_len_mpi, mpi_part, P_grid_right, 0, &(test.recv_buffer_west[0]), hor_len_mpi, mpi_part, P_grid_left, 0, grid_comm, &status);

    test.write_output_buffer(3, P_grid_rank);
    test.update_part();
    test.write_output_vec(3, P_grid_rank);

    MPI_Finalize();

    return 0;
}
