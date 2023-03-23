#include "species.hh"
#include "math.h"

species::species(std::string name_a, int *ppc_a, int *range_a, float *vf_a, float *vth_a) : name(name_a)
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

    N_x = range[0] + 1;
    N_y = range[1] + 1;

    xbox = dx * N_x;
    ybox = dy * N_y;

    // reserve space for the arrays of particles
    part A;
    vec.reserve(3 * np);
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

    // charge field initialization
    charge = new FCPIC::field();

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

    std::cout << "ppc0: " << ppc[0] << "  ppc[1]: " << ppc[1] << std::endl;

    loccell.reserve(np);

    for (int j = 0; j < ppc[1]; j++)
    {
        for (int i = 0; i < ppc[0]; i++)
        {
            loccell.push_back(dpcellx * ((float)i + 0.5f)); // In the middle of each subdivision
            loccell.push_back(dpcelly * ((float)j + 0.5f));
        }
    }

    // debugging
    // for (int i = 0; i < loccell.size(); i++)
    //     cout << "i: " << i << " value: " << loccell[i] << endl;

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
    std::vector<double> charge_vec;
    // without ghost cells
    charge_vec.assign(((N_x) * (N_y)), 0.f);
    for (int i = 0; i < np; i++)
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
    charge->val = charge_vec;
    charge->N_x = N_x;
    charge->N_y = N_y;
    charge->N = N_x * N_y;

    charge_vec.clear();
}

// dummy method to simulate the change of particles
// to be integrated in the pusher method
// void species::to_buffer()
// {
//     send_buffer_north.clear();
//     send_buffer_south.clear();
//     send_buffer_east.clear();
//     send_buffer_west.clear();

//     //! just an assumptoion to test....
//     int north_aux = (int)np / 2;
//     int south_aux = (int)np / 3;
//     int east_aux = (int)np / 4;
//     int west_aux = (int)np / 5;

//     for (int i = 0; i < buffer_north_len; i++)
//     {
//         int n_aux = north_aux + i;
//         vec[n_aux].x = -2;
//         send_buffer_north.push_back(vec[n_aux]);
//         // std::cout << vec[n_aux].ix << "  " << vec[n_aux].iy << std::endl;
//     }

//     for (int i = 0; i < buffer_south_len; i++)
//     {
//         int s_aux = south_aux + i;
//         vec[s_aux].x = -2;
//         send_buffer_south.push_back(vec[s_aux]);
//         // std::cout << vec[s_aux].ix << "  " << vec[s_aux].iy << std::endl;
//     }

//     for (int i = 0; i < buffer_east_len; i++)
//     {
//         int e_aux = east_aux + i;
//         vec[e_aux].x = -2;
//         send_buffer_east.push_back(vec[e_aux]);
//         // std::cout << vec[e_aux].ix << "  " << vec[e_aux].iy << std::endl;
//     }
//     for (int i = 0; i < buffer_west_len; i++)
//     {
//         int w_aux = west_aux + i;
//         vec[w_aux].x = -2;
//         send_buffer_west.push_back(vec[w_aux]);
//         // std::cout << vec[w_aux].ix << "  " << vec[w_aux].iy << std::endl;
//     }
//     //! just an assumption to test....

//     // erase elements out of box dimension
//     vec.erase(std::remove_if(vec.begin(), vec.end(), [this](const part o)
//                              { return (o.x < 0. || o.x > this->xbox || o.y < 0. || o.y > this->ybox); }),
//               vec.end());

//     //     // // debuggung buffers
//     //     // for (int i = 0; i < vec_buffer_north.size(); i++)
//     //     // {
//     //     //     std::cout << vec_buffer_north[i].ix << " " << vec_buffer_north[i].iy << " " << vec_buffer_north[i].x << "  " << vec_buffer_north[i].y << std::endl;
//     //     // }

//     //     //! Set the recv_buffer ix and iy to -1 to know whether or not there was exchange information in the MPI
//     part recv_dummy;
//     recv_dummy.ix = -1;
//     recv_dummy.iy = -1;
//     // actual memory allocation before passing messages through MPI
//     recv_buffer_north.assign(buffer_south_len, recv_dummy); //! assumption; improved later
//     recv_buffer_south.assign(buffer_north_len, recv_dummy);
//     recv_buffer_east.assign(buffer_west_len, recv_dummy);
//     recv_buffer_west.assign(buffer_east_len, recv_dummy);
// }

// void species::update_part()
// {
//     for (int i = 0; i < buffer_south_len; i++)
//     {
//         if (recv_buffer_north[i].ix != -1) // checking if there was "meaningful" comunication
//             vec.push_back(recv_buffer_north[i]);
//     }
//     for (int i = 0; i < buffer_south_len; i++)
//     {
//         if (recv_buffer_south[i].ix != -1)
//             vec.push_back(recv_buffer_south[i]);
//     }
//     for (int i = 0; i < buffer_west_len; i++)
//     {
//         if (recv_buffer_east[i].ix != -1)
//             vec.push_back(recv_buffer_east[i]);
//     }
//     for (int i = 0; i < buffer_east_len; i++)
//     {
//         if (recv_buffer_west[i].ix != -1)
//             vec.push_back(recv_buffer_west[i]);
//     }

//     recv_buffer_north.clear();
//     recv_buffer_south.clear();
//     recv_buffer_east.clear();
//     recv_buffer_west.clear();
// }

//???????????????????????????? later
void species::init_pusher(const float Ex, const float Ey, const int counter)
{
    int N_part = vec.size();

    // for (int counter = 0; counter < N_part; counter++)
    // {
    // float Ex_field = 0.f;
    // float Ey_field = 0.f;

    // get E field in the particle position
    //  E_field->get_field(vec[counter], Ex_field, Ey_field);
    float dt = 1.;

    vec[counter].ux = vec[counter].ux - 0.5 * q / m * Ex * dt;
    vec[counter].uy = vec[counter].uy - 0.5 * q / m * Ey * dt;
    // }
}

void species::particle_pusher(const float Ex, const float Ey, const int counter)
{
    int N_part = vec.size();

    // for (int counter = 0; counter < N_part; counter++)
    // {
    //     float Ex_field = 0.f;
    //     float Ey_field = 0.f;

    // get E field in the particle position
    //  E_field->get_field(vec[counter], Ex_field, Ey_field);
    float dt = 1.;

    vec[counter].ux = vec[counter].ux + q / m * Ex * dt;
    vec[counter].uy = vec[counter].uy + q / m * Ey * dt;

    vec[counter].x = vec[counter].x + vec[counter].ux * dt;
    vec[counter].y = vec[counter].y + vec[counter].uy * dt;

    //     // check boundaries cell in species grid
    // }
}

void species::advance_cell(int counter, int flag_MPI)
{
    float dx = 1.f;
    float dy = 1.f;

    float posx = vec[counter].x;
    float posy = vec[counter].y;

    bool xmin_cond = posx < 0.f;
    bool xmax_cond = posx > dx;
    bool ymin_cond = posy < 0.f;
    bool ymax_cond = posy > dy;

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
            // std::cout << "here" << std::endl;
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

    std::cout << "mid advance ****************" << std::endl;
    std::cout << "cell (" << vec[counter].ix << ", " << vec[counter].iy << ")" << std::endl;
    std::cout << "x:" << vec[counter].x << std::endl;
    std::cout << "y: " << vec[counter].y << std::endl;
    std::cout << "ux: " << vec[counter].ux << std::endl;
    std::cout << "uy: " << vec[counter].uy << std::endl;
    std::cout << "---------- ****************" << std::endl;

    // check if particle is going to be passed for MPI
    bool ixmin_cond = vec[counter].ix <= -1;
    bool ixmax_cond = vec[counter].ix >= range[0];
    bool iymin_cond = vec[counter].iy <= -1;
    bool iymax_cond = vec[counter].iy >= range[1];

    if (flag_MPI != -2) // periodic
    {
        while (vec[counter].ix < 0)
            vec[counter].ix = vec[counter].ix + range[0];

        while (vec[counter].iy < 0)
            vec[counter].iy = vec[counter].iy + range[1];

        while (vec[counter].ix >= range[0])
            vec[counter].ix = vec[counter].ix - range[0];

        while (vec[counter].iy >= range[1])
            vec[counter].iy = vec[counter].iy - range[1];

        if ((!ixmin_cond) && (!ixmax_cond) && iymax_cond)
        {
            std::cout << "north" << std::endl;
            send_buffer_north.push_back(vec[counter]);
        }
        else if (ixmax_cond && iymax_cond)
        {
            std::cout << "ne" << std::endl;
            send_buffer_ne.push_back(vec[counter]);
        }
        else if (ixmin_cond && iymax_cond)
        {
            std::cout << "nw" << std::endl;
            send_buffer_nw.push_back(vec[counter]);
        }
        // south buffer
        else if ((!ixmin_cond) && (!ixmax_cond) && iymin_cond)
        {
            std::cout << "south" << std::endl;
            send_buffer_south.push_back(vec[counter]);
        }
        else if (ixmax_cond && iymin_cond)
        {
            std::cout << "se" << std::endl;
            send_buffer_se.push_back(vec[counter]);
        }
        else if (ixmin_cond && iymin_cond)
        {
            std::cout << "sw" << std::endl;
            send_buffer_sw.push_back(vec[counter]);
        }
        // east buffer (vec[counter].ix >= range[1] && vec[counter].iy > -1)
        else if (ixmax_cond && (!iymin_cond) && (!iymax_cond))
        {
            std::cout << "east" << std::endl;
            send_buffer_east.push_back(vec[counter]);
        }
        // west buffer vec[counter].ix <= -1 && vec[counter].iy > -1
        else if (ixmin_cond && (!iymin_cond) && (!iymax_cond))
        {
            std::cout << "west" << std::endl;
            send_buffer_west.push_back(vec[counter]);
        }
        else
        {
            std::cout << "Nothing happens" << std::endl;
        }
    }

    if (flag_MPI == -2) // non periodic - physical
    {
        while (vec[counter].ix < -range[0])
            vec[counter].ix = vec[counter].ix + range[0];

        while (vec[counter].iy < -range[1])
            vec[counter].iy = vec[counter].iy + range[1];

        while (vec[counter].ix > range[0])
            vec[counter].ix = vec[counter].ix - range[0];

        while (vec[counter].iy > range[1])
            vec[counter].iy = vec[counter].iy - range[1];

        vec[counter].ix = fabs(vec[counter].ix);
        vec[counter].iy = fabs(vec[counter].iy);
    }

    //!!! change considitions to send
    // north buffer
    // if ((!ixmin_cond) && imax_cond)
    // {
    //     std::cout << "north" << std::endl;
    //     if (flag_MPI != -2)
    //     {
    //         send_buffer_north.push_back(vec[counter]);
    //     }
    //     else // physical boundaries
    //     {
    //         vec[counter].iy = range[1] - 1;
    //         vec[counter].y = dy;
    //         vec[counter].uy = -vec[counter].uy;
    //     }
    // }
    // // south buffer
    // // (!ixmin_cond) && ()vec[counter].iy <= -1
    // else if ((!ixmin_cond) && iymin_cond)
    // {
    //     std::cout << "south" << std::endl;
    //     if (flag_MPI == -2)
    //         send_buffer_south.push_back(vec[counter]);
    //     else // physical boundaries
    //     {
    //         vec[counter].iy = range[0];
    //         vec[counter].y = 0.;
    //         vec[counter].uy = -vec[counter].uy;
    //     }
    // } // east buffer (vec[counter].ix >= range[1] && vec[counter].iy > -1)
    // else if (imax_cond && (!iymin_cond))
    // {
    //     std::cout << "east" << std::endl;
    //     if (flag_MPI == -2)
    //         send_buffer_east.push_back(vec[counter]);
    //     else
    //     {
    //         vec[counter].ix = range[1] - 1;
    //         vec[counter].x = dx;
    //         vec[counter].ux = -vec[counter].ux;
    //     }
    // }
    // // west buffer vec[counter].ix <= -1 && vec[counter].iy > -1
    // else if (ixmin_cond && (!iymin_cond))
    // {
    //     std::cout << "west" << std::endl;
    //     if (flag_MPI == -2)
    //         send_buffer_west.push_back(vec[counter]);
    //     else
    //     {
    //         vec[counter].ix = range[0];
    //         vec[counter].x = 0.;
    //         vec[counter].ux = -vec[counter].ux;
    //     }
    // }
    // else
    // {
    //     std::cout << "Nothing happens" << std::endl;
    // }
}
// after loop over all particles

void species::to_buffer()
{
    // clean all particles
    vec.erase(std::remove_if(vec.begin(), vec.end(), [this](const part obj)
                             { return (obj.ix < 0 || obj.ix >= this->range[1] || obj.iy < 0 || obj.iy >= this->range[1]); }),
              vec.end());

    part recv_dummy;
    recv_dummy.ix = -1;
    recv_dummy.iy = -1;
    // actual memory allocation before passing messages through MPI
    recv_buffer_north.assign(buffer_south_len, recv_dummy); //! assumption; improved later
    recv_buffer_south.assign(buffer_north_len, recv_dummy);
    recv_buffer_east.assign(buffer_west_len, recv_dummy);
    recv_buffer_west.assign(buffer_east_len, recv_dummy);
}

// after communication
void species::update_part()
{
    for (int i = 0; i < buffer_south_len; i++)
    {
        if (recv_buffer_north[i].ix != -1) // checking if there was "meaningful" comunication
            vec.push_back(recv_buffer_north[i]);
    }
    for (int i = 0; i < buffer_south_len; i++)
    {
        if (recv_buffer_south[i].ix != -1)
            vec.push_back(recv_buffer_south[i]);
    }
    for (int i = 0; i < buffer_west_len; i++)
    {
        if (recv_buffer_east[i].ix != -1)
            vec.push_back(recv_buffer_east[i]);
    }
    for (int i = 0; i < buffer_east_len; i++)
    {
        if (recv_buffer_west[i].ix != -1)
            vec.push_back(recv_buffer_west[i]);
    }

    // clear old data from buffer's
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

//???????????????????????????? later

// for debugging
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

void species::write_output_vec(int rank, int time)
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
        Output_file << std::endl;
    }
    Output_file << std::endl
                << std::endl;
    for (int i = 0; i < charge->val.size(); i++)
    {
        if (i % (N_x) == 0)
            Output_file << std::endl;

        Output_file << std::setw(precision) << charge->val[i] << space;
        // if (i % 4 == 0)
    }

    Output_file.close();
}

void species::write_output_buffer(int rank, int time)
{
    std::fstream Output_file;
    std::string filename;

    filename = "../results/buffer_" + std::to_string(time) + "__t_" + std::to_string(time) + ".txt";
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
        Output_file << std::endl;
    }

    Output_file << std::endl
                << std::endl;

    Output_file.close();
}