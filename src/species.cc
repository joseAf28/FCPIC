#include "species.hh"
#include "math.h"
#include "mpi.h"

species::species(std::string name_a, int *ppc_a, int *range_a, float *vf_a, float *vth_a, float charge) : name(name_a), q(charge)
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

    // std::cout << "vth[0]: " << vth[0] << " vth[1]: " << vth[1] << " vth[2]: " << vth[2] << std::endl;

    // number of cells in each direction of the process domain
    N_int_x = range[0] + 1;
    N_int_y = range[1] + 1;

    N_x = range[0] + 3;
    N_y = range[1] + 3;

    // initializing vector with set_np_part() number: number of particles
    np = (N_int_x + 1) * ppc[0] * (N_int_y + 1) * ppc[1];

    // reserve space for the arrays of particles
    part A;
    A.flag = BULK;
    vec.reserve(3 * np); // assumption for the minimum reserved space
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

    // std::cout << __PRETTY_FUNCTION__ << std::endl;
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

    // std::cout << __PRETTY_FUNCTION__ << std::endl;
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
    for (int j = 0; j <= N_int_x; j++)
    {
        for (int i = 0; i <= N_int_y; i++)
        {
            for (int k = 0; k < npcell; k++)
            {
                vec[ip].ix = j;
                vec[ip].iy = i;
                vec[ip].x = loccell[2 * k];
                vec[ip].y = loccell[2 * k + 1];
                ip = ip + 1;
            }
        }
    }
    loccell.clear();
}

void species::get_charge(FCPIC::field *charge)
{
    //! Done in simulation: avoid creating a new field and doing the sum of all components
    // charge->setValue(0.f);

    int i, j;
    float wx, wy;
    for (int k = 0; k < vec.size(); k++)
    {
        i = vec[k].iy;
        j = vec[k].ix;
        wx = vec[k].x;
        wy = vec[k].y;

        // if (i > range[1])
        //     i = range[1];

        // if (j > range[0])
        //     j = range[0];

        // if (i < 0)
        //     i = 0;

        // if (j < 0)
        //     j = 0;

        charge->val[POSITION] += (dx - wx) * (dy - wy) * q / (dx * dy);
        charge->val[EAST] += wx * (dy - wy) * q / (dx * dy);
        charge->val[NORTH] += (dx - wx) * wy * q / (dx * dy);
        charge->val[NORTHEAST] += wx * wy * q / (dx * dy);
    }

    // TO BE UPDATED WITH N0
}

void species::field_inter(FCPIC::field *Ex, FCPIC::field *Ey, float &Ex_i, float &Ey_i, int counter)
{
    int i = vec[counter].iy;
    int j = vec[counter].ix;

    float wx = vec[counter].x;
    float wy = vec[counter].y;

    float A_pos = (dx - wx) * (dy - wy);
    float A_e = wx * (dy - wy);
    float A_n = (dx - wx) * wy;
    float A_ne = wx * wy;

    // if (i > range[1])
    //     i = range[1];

    // if (j > range[0])
    //     j = range[0];

    // if (i < 0)
    //     i = 0;

    // if (j < 0)
    //     j = 0;

    Ex_i = A_pos * Ex->val[POSITION] + A_e * Ex->val[EAST] + A_n * Ex->val[NORTH] + A_ne * Ex->val[NORTHEAST];

    Ey_i = A_pos * Ey->val[POSITION] + A_e * Ey->val[EAST] + A_n * Ey->val[NORTH] + A_ne * Ey->val[NORTHEAST];

    Ex_i /= (dx * dy);
    Ey_i /= (dx * dy);
}

void species::init_pusher(FCPIC::field *Ex, FCPIC::field *Ey)
{

    for (int i = 0; i < np; i++)
    {
        float Ex_i = 0.f;
        float Ey_i = 0.f;

        field_inter(Ex, Ey, Ex_i, Ey_i, i);

        vec[i].ux = vec[i].ux - 0.5 * q / m * Ex_i * dt;
        vec[i].uy = vec[i].uy - 0.5 * q / m * Ey_i * dt;
    }
}

void species::particle_pusher(FCPIC::field *Ex, FCPIC::field *Ey)
{

    for (int i = 0; i < np; i++)
    {
        float Ex_i = 0.f;
        float Ey_i = 0.f;

        field_inter(Ex, Ey, Ex_i, Ey_i, i);

        vec[i].ux = vec[i].ux + q / m * Ex_i * dt;
        vec[i].uy = vec[i].uy + q / m * Ey_i * dt;

        vec[i].x = vec[i].x + vec[i].ux * dt;
        vec[i].y = vec[i].y + vec[i].uy * dt;
    }
}

bool species::advance_cell(int *ranks_mpi)
{ // ranks_mpi[0] - rank, ranks_mpi[1] - top, ranks_mpi[2] - bottom,
    //  ranks_mpi[3] - right, ranks_mpi[4] - left
    bool changes_made = false;
    for (int counter = 0; counter < np; counter++)
    {
        if (vec[counter].x >= 0 && vec[counter].x < dx && vec[counter].y >= 0 && vec[counter].y < dx && vec[counter].ix >= 0 && vec[counter].ix < range[0] && vec[counter].iy >= 0 && vec[counter].iy < range[1])
            continue;

        changes_made = true;

        vec[counter].ix += floor(vec[counter].x / dx);
        vec[counter].x = fmod(vec[counter].x, dx);

        if (vec[counter].x < 0)
            vec[counter].x += dx;

        vec[counter].iy += floor(vec[counter].y / dy);
        vec[counter].y = fmod(vec[counter].y, dy);

        if (vec[counter].y < 0)
            vec[counter].y += dy;

        bool flag_top = ranks_mpi[1] == MPI_PROC_NULL;
        bool flag_bottom = ranks_mpi[2] == MPI_PROC_NULL;
        bool flag_right = ranks_mpi[3] == MPI_PROC_NULL;
        bool flag_left = ranks_mpi[4] == MPI_PROC_NULL;

        bool ixmin_cond = vec[counter].ix <= -1;
        bool ixmax_cond = vec[counter].ix > range[0];
        bool iymin_cond = vec[counter].iy <= -1;
        bool iymax_cond = vec[counter].iy > range[1];

        bool send_N_true = false;
        bool send_S_true = false;
        bool send_W_true = false;
        bool send_E_true = false;

        if (ixmin_cond)
        {
            if (flag_left)
            {
                vec[counter].x = dx - vec[counter].x;
                // vec[counter].ix = -vec[counter].ix - 1;
                while (vec[counter].ix < 0)
                    vec[counter].ix = vec[counter].ix + range[0];

                vec[counter].ix = range[0] - vec[counter].ix;

                vec[counter].ix = 0;
                vec[counter].ux = fabs(vec[counter].ux);

                if (vec[counter].ix < 0)
                {
                    std::cout << "flag_left ix: " << vec[counter].ix << std::endl;
                }
            }
            else
            {
                while (vec[counter].ix < 0)
                    vec[counter].ix = vec[counter].ix + range[0];
                // vec[counter].ix = fabs(vec[counter].ix + N_int_x) - 1;
                send_W_true = true;
                if (vec[counter].ix < 0)
                {
                    std::cout << "flag_left guard ix: " << vec[counter].ix << std::endl;
                }
            }
        }

        if (iymin_cond)
        {
            if (flag_bottom)
            {
                vec[counter].y = dy - vec[counter].y;
                // vec[counter].iy = -vec[counter].iy - 1;

                while (vec[counter].iy < 0)
                    vec[counter].iy = vec[counter].iy + range[1];

                vec[counter].iy = range[1] - vec[counter].iy;

                vec[counter].uy = fabs(vec[counter].uy);

                if (vec[counter].iy < 0)
                {
                    std::cout << "flag_bottom iy: " << vec[counter].iy << std::endl;
                }
            }
            else
            {
                while (vec[counter].iy < 0)
                    vec[counter].iy = vec[counter].iy + range[1];

                // vec[counter].iy = fabs(vec[counter].iy + N_int_y) - 1;
                send_S_true = true;

                if (vec[counter].iy < 0)
                {
                    std::cout << "flag_bottom guard iy: " << vec[counter].iy << std::endl;
                }
            }
        }

        if (ixmax_cond)
        {
            if (flag_right)
            {
                vec[counter].x = dx - vec[counter].x;
                // vec[counter].ix = 2 * N_int_x - vec[counter].ix;

                while (vec[counter].ix > range[0])
                    vec[counter].ix = vec[counter].ix - range[0];

                vec[counter].ix = range[0] - vec[counter].ix;

                vec[counter].ux = -fabs(vec[counter].ux);

                if (vec[counter].ix < 0)
                {
                    std::cout << "flag_right ix: " << vec[counter].ix << std::endl;
                }
            }
            else
            {
                while (vec[counter].ix > range[0])
                    vec[counter].ix = vec[counter].ix - range[0];

                // vec[counter].ix = fabs(vec[counter].ix - N_int_x) + 1;
                send_E_true = true;
                if (vec[counter].ix < 0)
                {
                    std::cout << "flag_right guard ix: " << vec[counter].ix << std::endl;
                }
            }
        }

        if (iymax_cond)
        {
            if (flag_top)
            {
                vec[counter].y = dy - vec[counter].y;
                // vec[counter].iy = 2 * N_int_y - vec[counter].iy;

                while (vec[counter].iy > range[1])
                    vec[counter].iy = vec[counter].iy - range[1];

                vec[counter].iy = range[1] - vec[counter].iy;

                vec[counter].uy = -fabs(vec[counter].uy);

                if (vec[counter].iy < 0)
                {
                    std::cout << "flag_top iy: " << vec[counter].iy << std::endl;
                }
            }
            else
            {
                while (vec[counter].iy > range[1])
                    vec[counter].iy = vec[counter].iy - range[1];

                // vec[counter].iy = fabs(vec[counter].iy - N_int_y) + 1;
                send_N_true = true;

                if (vec[counter].iy < 0)
                {
                    std::cout << "flag_top guard iy: " << vec[counter].iy << std::endl;
                }
            }
        }

        // Periodic and virtual boundary conditions
        // north buffer
        if ((!send_W_true) && (!send_E_true) && send_N_true)
        {
            send_buffer_north.push_back(vec[counter]);
            vec[counter].flag = SEND;
        }
        // ne buffer
        else if (send_N_true && send_E_true)
        {
            send_buffer_ne.push_back(vec[counter]);
            vec[counter].flag = SEND;
        }
        // nw buffer
        else if (send_N_true && send_W_true)
        {
            send_buffer_nw.push_back(vec[counter]);
            vec[counter].flag = SEND;
        }
        // south buffer
        else if ((!send_W_true) && (!send_E_true) && send_S_true)
        {
            send_buffer_south.push_back(vec[counter]);
            vec[counter].flag = SEND;
        }
        // se buffer
        else if (send_S_true && send_E_true)
        {
            send_buffer_se.push_back(vec[counter]);
            vec[counter].flag = SEND;
        }
        // sw buffer
        else if (send_S_true && send_W_true)
        {
            send_buffer_sw.push_back(vec[counter]);
            vec[counter].flag = SEND;
        }
        // east buffer
        else if (send_E_true && (!send_N_true) && (!send_S_true))
        {
            send_buffer_east.push_back(vec[counter]);
            vec[counter].flag = SEND;
        }
        // west
        else if (send_W_true && (!send_N_true) && (!send_S_true))
        {
            send_buffer_west.push_back(vec[counter]);
            vec[counter].flag = SEND;
        }
        else
        {
        }
    }
    return changes_made;
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
    //
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

void species::write_output_vec(const int rank, const int spec, const int counter)
{
    std::fstream Output_file;
    std::string filename;

    filename = "../results/phase_space/particles_" + std::to_string(rank) + "_spec_" + std::to_string(spec) + "_counter_" + std::to_string(counter) + ".txt";
    std::string space = "   ";

    Output_file.open(filename, std::ios::out);
    Output_file << "x:  y:  ux:  uy:" << std::endl;

    int precision = 4;
    for (int i = 0; i < vec.size(); i++)
    {
        float posx = vec[i].ix * dx + vec[i].x;
        float posy = vec[i].iy * dy + vec[i].y;

        Output_file << posx << space;
        Output_file << posy << space;
        Output_file << vec[i].ux << space << vec[i].uy << std::endl;

        // Output_file << "cell (" << vec[i].ix << ", " << vec[i].iy << ")" << space;
        // Output_file << "x:" << std::setw(precision) << vec[i].x << space;
        // Output_file << "y: " << std::setw(precision) << vec[i].y << space;
        // Output_file << "ux: " << std::setw(precision) << vec[i].ux << space;
        // Output_file << "uy: " << std::setw(precision) << vec[i].uy << space;
        // Output_file << "flag: " << std::setw(precision) << vec[i].flag << space;
        // Output_file << std::endl;
    }
    Output_file << std::endl
                << std::endl;
    /*
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
    */
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
