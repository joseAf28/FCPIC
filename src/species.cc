#include "species.hh"
#include "math.h"
#include "mpi.h"

species::species(std::string name_a, int *ppc_a, int *range_a, float *vf_a, float *vth_a) : name(name_a), charge(nullptr)
{
    ppc[0] = ppc_a[0];
    ppc[1] = ppc_a[1];

    range[0] = range_a[0];
    range[1] = range_a[1];

    vf[0] = vf_a[0];
    vf[1] = vf_a[1];
    vf[2] = vf_a[2];

    vth[0] = vth_a[0];
    vth[1] = vth_a[1];
    vth[2] = vth_a[2];

    // initializing vector with set_np_part() number: number of particles
    np = range[0] * ppc[0] * range[1] * ppc[1];

    // number of cells in each direction of the process domain
    N_x = range[0] + 1;
    N_y = range[1] + 1;

    // reserve space for the arrays of particles
    part A;
    A.flag = BULK;
    vec.reserve(3 * np); // assumotion for the minimum reserved space
    vec.assign(np, A);
    send_buffer_north.reserve(np); // assumption for the space
    send_buffer_south.reserve(np); // ! Think about it later
    send_buffer_east.reserve(np);
    send_buffer_west.reserve(np);

    send_buffer_ne.reserve(np); // assumption for the space
    send_buffer_se.reserve(np); // ! Think about it later
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

    // random number generator
    std::random_device dev;
    std::mt19937_64 rng_aux(dev());
    std::normal_distribution<double> norm_aux(0, 1);
    rng = rng_aux;
    rand_gauss = norm_aux;

    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

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

    delete charge;
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

void species::set_u()
{
    for (int i = 0; i < np; i++)
    {
        vec[i].ux = vf[0] + vth[0] * rand_gauss(rng);
        vec[i].uy = vf[1] + vth[1] * rand_gauss(rng);
        vec[i].uz = vf[2] + vth[2] * rand_gauss(rng);
    }
}

void species::set_x()
{
    std::vector<float> loccell;

    const int npcell = ppc[0] * ppc[1];
    const float dpcellx = dx / ppc[0];
    const float dpcelly = dy / ppc[1];

    loccell.reserve(np);

    for (int j = 0; j < ppc[1]; j++)
    {
        for (int i = 0; i < ppc[0]; i++)
        {
            loccell.push_back(dpcellx * ((float)i + 0.5f)); // In the middle of each subdivision
            loccell.push_back(dpcelly * ((float)j + 0.5f));
        }
    }

    int ip = 0;
    //! Uniform Density of Particles
    for (int j = 0; j < range[1]; j++)
    {
        for (int i = 0; i < range[0]; i++)
        {
            for (int k = 0; k < npcell; k++)
            {
                vec[ip].ix = i;
                vec[ip].iy = j;
                vec[ip].x = loccell[2 * k];
                vec[ip].y = loccell[2 * k + 1];
                ip = ip + 1;
            }
        }
    }
    loccell.clear();
}

void species::get_charge()
{
    //! charge vec with no ghost cells: ghost cells are going to be taken care in Fields class
    if (charge != nullptr)
        delete charge;

    std::vector<double> charge_vec((N_x) * (N_y), 0.f);

    for (int i = 0; i < vec.size(); i++)
    {
        int ij = vec[i].ix + N_x * vec[i].iy;

        float wx = vec[i].x;
        float wy = vec[i].y;

        // divide by dx*dy
        charge_vec[ij] += (dx - wx) * (dy - wy) * q;
        charge_vec[ij + 1] += wx * (dy - wy) * q;
        charge_vec[ij + N_x] += (dx - wx) * wy * q;
        charge_vec[ij + 1 + N_x] += wx * wy * q;
    }
    // update Field Charge
    charge = new FCPIC::field(N_x, N_y, charge_vec); // The physical domain, with no cell guards is sent

    charge_vec.clear();
}

//!! Define the Efield in each particle: missing feature
void species::init_pusher(FCPIC::field *Ex, FCPIC::field *Ey)
{
    int N_part = vec.size();

    for (int i = 0; i < vec.size(); i++)
    {
        //! Interpolation inside the cell is missing
        // get E field in the particle position
        int ij = vec[i].ix + N_x * vec[i].iy;

        vec[i].ux = vec[i].ux - 0.5 * q / m * Ex->val[ij] * dt;
        vec[i].uy = vec[i].uy - 0.5 * q / m * Ey->val[ij] * dt;
    }
}

//!! Define the Efield in each particle: missing feature
void species::particle_pusher(FCPIC::field *Ex, FCPIC::field *Ey)
{
    for (int i = 0; i < vec.size(); i++)
    {
        //! Interpolation inside the cell is missing
        // get E field in the particle position
        int ij = vec[i].ix + N_x * vec[i].iy;

        vec[i].ux = vec[i].ux + q / m * Ex->val[ij] * dt;
        vec[i].uy = vec[i].uy + q / m * Ey->val[ij] * dt;

        vec[i].x = vec[i].x + vec[i].ux * dt;
        vec[i].y = vec[i].y + vec[i].uy * dt;
    }
}

void species::advance_cell(int *ranks_mpi)
{ // ranks_mpi[0] - rank, ranks_mpi[1] - top, ranks_mpi[2] - bottom,
    //  ranks_mpi[3] - right, ranks_mpi[4] - left
    for (int counter = 0; counter < vec.size(); counter++)
    {
        float dx = 1.f;
        float dy = 1.f;

        float posx = vec[counter].x;
        float posy = vec[counter].y;

        bool xmin_cond = posx < 0.f;
        bool xmax_cond = posx > dx;
        bool ymin_cond = posy < 0.f;
        bool ymax_cond = posy > dy;

        // // Debugging motion
        // std::cout << "Init ****************" << std::endl;
        // std::cout << "cell (" << vec[counter].ix << ", " << vec[counter].iy << ")" << std::endl;
        // std::cout << "x:" << vec[counter].x << std::endl;
        // std::cout << "y: " << vec[counter].y << std::endl;
        // std::cout << "ux: " << vec[counter].ux << std::endl;
        // std::cout << "uy: " << vec[counter].uy << std::endl;
        // std::cout << "---------- ****************" << std::endl;

        while ((xmin_cond || xmax_cond || ymin_cond || ymax_cond))
        {
            if (xmin_cond)
            {
                vec[counter].x = dx - fabs(posx);
                vec[counter].ix = vec[counter].ix - 1;
            }
            if (xmax_cond)
            {
                vec[counter].x = posx - dx;
                vec[counter].ix = vec[counter].ix + 1;
            }
            if (ymin_cond)
            {
                vec[counter].y = dy - fabs(posy);
                vec[counter].iy = vec[counter].iy - 1;
            }
            if (ymax_cond)
            {
                vec[counter].y = posy - dy;
                vec[counter].iy = vec[counter].iy + 1;
            }
            posx = vec[counter].x;
            posy = vec[counter].y;

            xmin_cond = posx < 0.f;
            xmax_cond = posx > dx;
            ymin_cond = posy < 0.f;
            ymax_cond = posy > dy;
        }

        // // debugging the particles' motions inside the cell
        // std::cout << "mid  ****************" << std::endl;
        // std::cout << "cell (" << vec[counter].ix << ", " << vec[counter].iy << ")" << std::endl;
        // std::cout << "x:" << vec[counter].x << std::endl;
        // std::cout << "y: " << vec[counter].y << std::endl;
        // std::cout << "ux: " << vec[counter].ux << std::endl;
        // std::cout << "uy: " << vec[counter].uy << std::endl;
        // std::cout << "---------- ****************" << std::endl;

        bool flag_top = ranks_mpi[1] == MPI_PROC_NULL;
        bool flag_bottom = ranks_mpi[2] == MPI_PROC_NULL;
        bool flag_right = ranks_mpi[3] == MPI_PROC_NULL;
        bool flag_left = ranks_mpi[4] == MPI_PROC_NULL;

        int reflection = false; // check if we have reflection or not

        if (flag_top || flag_bottom || flag_right || flag_left) // non periodic - physical - boundaries
        {
            // std::cout << "grid_rank: " << ranks_mpi[0] << " reflection" << std::endl;
            // reflection of the particles
            if (flag_left && vec[counter].ix < 0)
            {
                vec[counter].ux = fabs(vec[counter].ux);
                vec[counter].ix = 0;
                vec[counter].x = 0.;
                vec[counter].flag = BULK;
                reflection = true;
            }

            if (flag_right && vec[counter].ix >= range[0])
            {
                vec[counter].ux = -fabs(vec[counter].ux);
                vec[counter].ix = range[0] - 1;
                vec[counter].x = dx;
                vec[counter].flag = BULK;
                reflection = true;
            }

            if (flag_bottom && vec[counter].iy < 0)
            {
                vec[counter].uy = fabs(vec[counter].uy);
                vec[counter].iy = 0;
                vec[counter].y = 0.;
                vec[counter].flag = BULK;
                reflection = true;
            }

            if (flag_top && vec[counter].iy >= range[1])
            {
                vec[counter].uy = -fabs(vec[counter].uy);
                vec[counter].iy = range[1] - 1;
                vec[counter].y = dx;
                vec[counter].flag = BULK;
                reflection = true;
            }

            // assunption: reflection takes priority over exchange particles with MPI
            if (vec[counter].ix < 0 && reflection == true)
                vec[counter].ix = 0;

            if (vec[counter].iy < 0 && reflection == true)
                vec[counter].iy = 0;

            if (vec[counter].ix >= range[0] && reflection == true)
                vec[counter].ix = range[0] - 1;

            if (vec[counter].iy >= range[1] && reflection == true)
                vec[counter].iy = vec[counter].iy - range[1];

            if (reflection == true)
                continue;
        }

        // std::cout << "end ****************" << std::endl;
        // std::cout << "cell (" << vec[counter].ix << ", " << vec[counter].iy << ")" << std::endl;
        // std::cout << "x:" << vec[counter].x << std::endl;
        // std::cout << "y: " << vec[counter].y << std::endl;
        // std::cout << "ux: " << vec[counter].ux << std::endl;
        // std::cout << "uy: " << vec[counter].uy << std::endl;
        // std::cout << "---------- ****************" << std::endl;

        // check if particle is going to be passed for MPI
        bool ixmin_cond = vec[counter].ix <= -1;
        bool ixmax_cond = vec[counter].ix >= range[0];
        bool iymin_cond = vec[counter].iy <= -1;
        bool iymax_cond = vec[counter].iy >= range[1];

        if (ranks_mpi[0] != MPI_PROC_NULL) // periodic
        {
            while (vec[counter].ix < 0)
                vec[counter].ix = vec[counter].ix + range[0];

            while (vec[counter].iy < 0)
                vec[counter].iy = vec[counter].iy + range[1];

            while (vec[counter].ix >= range[0])
                vec[counter].ix = vec[counter].ix - range[0];

            while (vec[counter].iy >= range[1])
                vec[counter].iy = vec[counter].iy - range[1];

            // north buffer
            if ((!ixmin_cond) && (!ixmax_cond) && iymax_cond)
            {
                // std::cout << "north" << std::endl;
                send_buffer_north.push_back(vec[counter]);
                vec[counter].flag = SEND; //!!!!mark to delete in prepare_to_buffer method
            }                             // ne buffer
            else if (ixmax_cond && iymax_cond)
            {
                // std::cout << "ne" << std::endl;
                send_buffer_ne.push_back(vec[counter]);
                vec[counter].flag = SEND;
            } // nw buffer
            else if (ixmin_cond && iymax_cond)
            {
                // std::cout << "nw" << std::endl;
                send_buffer_nw.push_back(vec[counter]);
                vec[counter].flag = SEND;
            } // south buffer
            else if ((!ixmin_cond) && (!ixmax_cond) && iymin_cond)
            {
                // std::cout << "south" << std::endl;
                send_buffer_south.push_back(vec[counter]);
                vec[counter].flag = SEND;
            } // se buffer
            else if (ixmax_cond && iymin_cond)
            {
                // std::cout << "se" << std::endl;
                send_buffer_se.push_back(vec[counter]);
                vec[counter].flag = SEND;
            } // sw buffer
            else if (ixmin_cond && iymin_cond)
            {
                // std::cout << "sw" << std::endl;
                send_buffer_sw.push_back(vec[counter]);
                vec[counter].flag = SEND;
            } // east buffer
            else if (ixmax_cond && (!iymin_cond) && (!iymax_cond))
            {
                // std::cout << "east" << std::endl;
                send_buffer_east.push_back(vec[counter]);
                vec[counter].flag = SEND;
            } // west
            else if (ixmin_cond && (!iymin_cond) && (!iymax_cond))
            {
                // std::cout << "west" << std::endl;
                send_buffer_west.push_back(vec[counter]);
                vec[counter].flag = SEND;
            }
            else
            {
                // std::cout << "Nothing happens" << std::endl;
            }
        }
    }
}

void species::prepare_buffer()
{
    // determine the size of the arrays that Exchange particles in MPI
    size_send_north = send_buffer_north.size();
    size_send_south = send_buffer_south.size();
    size_send_east = send_buffer_east.size();
    size_send_west = send_buffer_west.size();

    size_send_ne = send_buffer_ne.size();
    size_send_se = send_buffer_se.size();
    size_send_nw = send_buffer_nw.size();
    size_send_sw = send_buffer_sw.size();

    // delete all particles that were sent to exchange array buffers
    vec.erase(std::remove_if(vec.begin(), vec.end(), [this](const part obj)
                             { return (obj.flag == SEND); }),
              vec.end());
}

void species::update_part_list()
{
    // include the new particles after the MPI data exchange
    for (int i = 0; i < recv_buffer_north.size(); i++)
    {
        if (recv_buffer_north[i].ix != -1) // checking if there was real MPI communication
            vec.push_back(recv_buffer_north[i]);
    }
    for (int i = 0; i < recv_buffer_south.size(); i++)
    {
        if (recv_buffer_south[i].ix != -1)
            vec.push_back(recv_buffer_south[i]);
    }
    for (int i = 0; i < recv_buffer_east.size(); i++)
    {
        if (recv_buffer_east[i].ix != -1)
            vec.push_back(recv_buffer_east[i]);
    }
    for (int i = 0; i < recv_buffer_west.size(); i++)
    {
        if (recv_buffer_west[i].ix != -1)
            vec.push_back(recv_buffer_west[i]);
    }
    for (int i = 0; i < recv_buffer_ne.size(); i++)
    {
        if (recv_buffer_ne[i].ix != -1)
            vec.push_back(recv_buffer_ne[i]);
    }
    for (int i = 0; i < recv_buffer_nw.size(); i++)
    {
        if (recv_buffer_nw[i].ix != -1)
            vec.push_back(recv_buffer_nw[i]);
    }
    for (int i = 0; i < recv_buffer_sw.size(); i++)
    {
        if (recv_buffer_sw[i].ix != -1)
            vec.push_back(recv_buffer_sw[i]);
    }
    for (int i = 0; i < recv_buffer_se.size(); i++)
    {
        if (recv_buffer_se[i].ix != -1)
            vec.push_back(recv_buffer_se[i]);
    }

    // update the particles number in the simulation after MPI exchange
    np = vec.size();

    // clean old data from buffer's
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

//!!!!!!!!!!!!!!!!!! temporary methods for debugging methods
void species::print()
{
    for (int i = 0; i < np; i++)
    {
        std::cout << "cell (" << vec[i].ix << ", " << vec[i].iy << ")" << std::endl;
        std::cout << "x:" << vec[i].x << std::endl;
        std::cout << "y: " << vec[i].y << std::endl;
        std::cout << "ux: " << vec[i].ux << std::endl;
        std::cout << "uy: " << vec[i].uy << std::endl;
        std::cout << "****************" << std::endl;
    }
}

void species::write_output_vec(const int rank, const int time)
{
    std::fstream Output_file;
    std::string filename;

    filename = "../results/grid_" + std::to_string(rank) + "__t_" + std::to_string(time) + ".txt";
    std::string space = "   ";

    Output_file.open(filename, std::ios::out);
    Output_file << "**************Grid**************" << std::endl;

    int precision = 4;
    for (int i = 0; i < vec.size(); i++)
    {
        Output_file << "cell (" << vec[i].ix << ", " << vec[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << vec[i].x << space;
        Output_file << "y: " << std::setw(precision) << vec[i].y << space;
        Output_file << "ux: " << std::setw(precision) << vec[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << vec[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << vec[i].flag << space;
        Output_file << std::endl;
    }
    Output_file << std::endl
                << std::endl;

    if (charge != nullptr)
    {
        for (int i = 0; i < charge->val.size(); i++)
        {
            if (i % (N_x) == 0)
                Output_file << std::endl;

            Output_file << std::setw(precision) << charge->val[i] << space;
            // if (i % 4 == 0)
        }
    }

    Output_file.close();
}

void species::write_output_buffer(const int rank, const int time)
{
    std::fstream Output_file;
    std::string filename;

    filename = "../results/buffer_output_" + std::to_string(rank) + "_:_" + std::to_string(time) + ".txt";
    std::string space = "   ";

    Output_file.open(filename, std::ios::out);

    int precision = 4;

    Output_file << "**************Buffer North**************" << std::endl;
    for (int i = 0; i < recv_buffer_north.size(); i++)
    {
        Output_file << "cell (" << recv_buffer_north[i].ix << ", " << recv_buffer_north[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << recv_buffer_north[i].x << space;
        Output_file << "y: " << std::setw(precision) << recv_buffer_north[i].y << space;
        Output_file << "ux: " << std::setw(precision) << recv_buffer_north[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << recv_buffer_north[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << recv_buffer_north[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Buffer South**************" << std::endl;
    for (int i = 0; i < recv_buffer_south.size(); i++)
    {
        Output_file << "cell (" << recv_buffer_south[i].ix << ", " << recv_buffer_south[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << recv_buffer_south[i].x << space;
        Output_file << "y: " << std::setw(precision) << recv_buffer_south[i].y << space;
        Output_file << "ux: " << std::setw(precision) << recv_buffer_south[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << recv_buffer_south[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << recv_buffer_south[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Buffer East**************" << std::endl;
    for (int i = 0; i < recv_buffer_east.size(); i++)
    {
        Output_file << "cell (" << recv_buffer_east[i].ix << ", " << recv_buffer_east[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << recv_buffer_east[i].x << space;
        Output_file << "y: " << std::setw(precision) << recv_buffer_east[i].y << space;
        Output_file << "ux: " << std::setw(precision) << recv_buffer_east[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << recv_buffer_east[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << recv_buffer_east[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Buffer West**************" << std::endl;
    for (int i = 0; i < recv_buffer_west.size(); i++)
    {
        Output_file << "cell (" << recv_buffer_west[i].ix << ", " << recv_buffer_west[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << recv_buffer_west[i].x << space;
        Output_file << "y: " << std::setw(precision) << recv_buffer_west[i].y << space;
        Output_file << "ux: " << std::setw(precision) << recv_buffer_west[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << recv_buffer_west[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << recv_buffer_west[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Buffer NW**************" << std::endl;
    for (int i = 0; i < recv_buffer_nw.size(); i++)
    {
        Output_file << "cell (" << recv_buffer_nw[i].ix << ", " << recv_buffer_nw[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << recv_buffer_nw[i].x << space;
        Output_file << "y: " << std::setw(precision) << recv_buffer_nw[i].y << space;
        Output_file << "ux: " << std::setw(precision) << recv_buffer_nw[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << recv_buffer_nw[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << recv_buffer_nw[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Buffer SW**************" << std::endl;
    for (int i = 0; i < recv_buffer_sw.size(); i++)
    {
        Output_file << "cell (" << recv_buffer_sw[i].ix << ", " << recv_buffer_sw[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << recv_buffer_sw[i].x << space;
        Output_file << "y: " << std::setw(precision) << recv_buffer_sw[i].y << space;
        Output_file << "ux: " << std::setw(precision) << recv_buffer_sw[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << recv_buffer_sw[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << recv_buffer_sw[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Buffer SE**************" << std::endl;
    for (int i = 0; i < recv_buffer_se.size(); i++)
    {
        Output_file << "cell (" << recv_buffer_se[i].ix << ", " << recv_buffer_se[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << recv_buffer_se[i].x << space;
        Output_file << "y: " << std::setw(precision) << recv_buffer_se[i].y << space;
        Output_file << "ux: " << std::setw(precision) << recv_buffer_se[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << recv_buffer_se[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << recv_buffer_se[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Buffer NE**************" << std::endl;
    for (int i = 0; i < recv_buffer_ne.size(); i++)
    {
        Output_file << "cell (" << recv_buffer_ne[i].ix << ", " << recv_buffer_ne[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << recv_buffer_ne[i].x << space;
        Output_file << "y: " << std::setw(precision) << recv_buffer_ne[i].y << space;
        Output_file << "ux: " << std::setw(precision) << recv_buffer_ne[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << recv_buffer_ne[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << recv_buffer_ne[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << std::endl
                << std::endl;

    Output_file.close();
}

void species::write_input_buffer(const int rank, const int time)
{
    std::fstream Output_file;
    std::string filename;

    filename = "../results/buffer_input_" + std::to_string(rank) + "_:_" + std::to_string(time) + ".txt";
    std::string space = "   ";

    Output_file.open(filename, std::ios::out);

    int precision = 4;

    Output_file << "**************Real Buffer North**************" << std::endl;
    for (int i = 0; i < send_buffer_north.size(); i++)
    {
        Output_file << "cell (" << send_buffer_north[i].ix << ", " << send_buffer_north[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << send_buffer_north[i].x << space;
        Output_file << "y: " << std::setw(precision) << send_buffer_north[i].y << space;
        Output_file << "ux: " << std::setw(precision) << send_buffer_north[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << send_buffer_north[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << send_buffer_north[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Real Buffer South**************" << std::endl;
    for (int i = 0; i < send_buffer_south.size(); i++)
    {
        Output_file << "cell (" << send_buffer_south[i].ix << ", " << send_buffer_south[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << send_buffer_south[i].x << space;
        Output_file << "y: " << std::setw(precision) << send_buffer_south[i].y << space;
        Output_file << "ux: " << std::setw(precision) << send_buffer_south[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << send_buffer_south[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << send_buffer_south[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Real Buffer East**************" << std::endl;
    for (int i = 0; i < send_buffer_east.size(); i++)
    {
        Output_file << "cell (" << send_buffer_east[i].ix << ", " << send_buffer_east[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << send_buffer_east[i].x << space;
        Output_file << "y: " << std::setw(precision) << send_buffer_east[i].y << space;
        Output_file << "ux: " << std::setw(precision) << send_buffer_east[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << send_buffer_east[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << send_buffer_east[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Real Buffer West**************" << std::endl;
    for (int i = 0; i < send_buffer_west.size(); i++)
    {
        Output_file << "cell (" << send_buffer_west[i].ix << ", " << send_buffer_west[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << send_buffer_west[i].x << space;
        Output_file << "y: " << std::setw(precision) << send_buffer_west[i].y << space;
        Output_file << "ux: " << std::setw(precision) << send_buffer_west[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << send_buffer_west[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << send_buffer_west[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Real Buffer NW**************" << std::endl;
    for (int i = 0; i < send_buffer_nw.size(); i++)
    {
        Output_file << "cell (" << send_buffer_nw[i].ix << ", " << send_buffer_nw[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << send_buffer_nw[i].x << space;
        Output_file << "y: " << std::setw(precision) << send_buffer_nw[i].y << space;
        Output_file << "ux: " << std::setw(precision) << send_buffer_nw[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << send_buffer_nw[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << send_buffer_nw[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Real Buffer SW**************" << std::endl;
    for (int i = 0; i < send_buffer_sw.size(); i++)
    {
        Output_file << "cell (" << send_buffer_sw[i].ix << ", " << send_buffer_sw[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << send_buffer_sw[i].x << space;
        Output_file << "y: " << std::setw(precision) << send_buffer_sw[i].y << space;
        Output_file << "ux: " << std::setw(precision) << send_buffer_sw[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << send_buffer_sw[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << send_buffer_sw[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Real Buffer SE**************" << std::endl;
    for (int i = 0; i < send_buffer_se.size(); i++)
    {
        Output_file << "cell (" << send_buffer_se[i].ix << ", " << send_buffer_se[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << send_buffer_se[i].x << space;
        Output_file << "y: " << std::setw(precision) << send_buffer_se[i].y << space;
        Output_file << "ux: " << std::setw(precision) << send_buffer_se[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << send_buffer_se[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << send_buffer_se[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << "**************Real Buffer NE**************" << std::endl;
    for (int i = 0; i < send_buffer_ne.size(); i++)
    {
        Output_file << "cell (" << send_buffer_ne[i].ix << ", " << send_buffer_ne[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << send_buffer_ne[i].x << space;
        Output_file << "y: " << std::setw(precision) << send_buffer_ne[i].y << space;
        Output_file << "ux: " << std::setw(precision) << send_buffer_ne[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << send_buffer_ne[i].uy << space;
        Output_file << "flag: " << std::setw(precision) << send_buffer_ne[i].flag << space;
        Output_file << std::endl;
    }

    Output_file << std::endl
                << std::endl;

    Output_file.close();
}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
