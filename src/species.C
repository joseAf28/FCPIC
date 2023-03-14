#include "species.h"

// declare all stuff
species::species(std::string name_a, int *ppc_a, int *range_a, int *vec_U_a, float *vf_a, float *vth_a)
{
    // declare all variables
    name = name_a;

    ppc[0] = ppc_a[0];
    ppc[1] = ppc_a[1];

    range[0] = range_a[0];
    range[1] = range_a[1];
    range[2] = range_a[2];
    range[3] = range_a[3];

    init_U = vec_U_a[0];
    end_U = vec_U_a[1];
    vf[0] = vf_a[0];
    vf[1] = vf_a[1];
    vf[2] = vf_a[2];

    vth[0] = vth_a[0];
    vth[1] = vth_a[1];
    vth[2] = vth_a[2];

    // initializing vector with set_np_part() number: number of particles
    np = set_nb();

    nx = (range[1] - range[0]) + 1;
    ny = (range[3] - range[2]) + 1;

    xbox = dx * nx;
    ybox = dy * ny;

    // initializing the array of particles
    vec = std::unique_ptr<part>(new part[np]);
    charge = new simulation::field(nx + 2, ny + 2);

    // random number generator
    std::random_device dev;
    std::mt19937_64 rng_aux(dev());
    std::normal_distribution<double> norm_aux(0, 1);
    rng = rng_aux;
    rand_gauss = norm_aux;

    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

species::~species() //! fix problem with deallocation of mat memory
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

int species::set_nb()
{
    // Assuming Uniform Density
    float nb_part = (range[1] - range[0]) * ppc[0] * (range[3] - range[2]) * ppc[1];
    return nb_part;
}

void species::set_U()
{
    end_U = np;
    for (int i = 0; i < init_U; i++)
    {
        vec.get()[i].ux = vf[0] + vth[0] * rand_gauss(rng);
        vec.get()[i].uy = vf[1] + vth[1] * rand_gauss(rng);
        vec.get()[i].uz = vf[2] + vth[2] * rand_gauss(rng);
    }
    for (int i = init_U; i < end_U; i++)
    {
        vec.get()[i].ux = vf[0] + vth[0] * rand_gauss(rng);
        vec.get()[i].uy = vf[1] + vth[1] * rand_gauss(rng);
        vec.get()[i].uz = vf[2] + vth[2] * rand_gauss(rng);
    }
    for (int i = end_U; i < np; i++)
    {
        vec.get()[i].ux = vf[0] + vth[0] * rand_gauss(rng);
        vec.get()[i].uy = vf[1] + vth[1] * rand_gauss(rng);
        vec.get()[i].uz = vf[2] + vth[2] * rand_gauss(rng);
    }
}

void species::set_X()
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

    // debugging
    // for (int i = 0; i < loccell.size(); i++)
    //     cout << "i: " << i << " value: " << loccell[i] << endl;

    int ip = 0;
    //! Uniform Density of Particles
    for (int j = range[2]; j < range[3]; j++)
    {
        for (int i = range[0]; i < range[1]; i++)
        {
            for (int k = 0; k < npcell; k++)
            {
                vec.get()[ip].ix = i;
                vec.get()[ip].iy = j;
                vec.get()[ip].x = loccell[2 * k];
                vec.get()[ip].y = loccell[2 * k + 1];
                ip = ip + 1;
            }
        }
    }
    loccell.clear();
}

void species::get_charge()
{
    std::vector<double> charge_vec(((nx + 2) * (ny + 2)), 0.f);

    // using matrices
    mat = new float *[ny + 2];
    for (int i = 0; i < nx + 2; i++)
        mat[i] = new float[nx + 2]();

    for (int i = 0; i < np; i++)
    {
        int ix = vec.get()[i].ix;
        int iy = vec.get()[i].iy;

        ix = ix + 1;
        iy = iy + 1;

        float wx = vec.get()[i].x;
        float wy = vec.get()[i].y;

        mat[ix][iy] += (dx - wx) * (dy - wy) * q;
        mat[ix + 1][iy] += wx * (dy - wy) * q;
        mat[ix][iy + 1] += (dx - wx) * wy * q;
        mat[ix + 1][iy + 1] += wx * wy * q;
    }
    // fill guard cells
    for (int j = 0; j < ny + 2; j++)
    {
        mat[0][j] = mat[1][j];
        mat[nx + 1][j] = mat[nx][j];
    }
    for (int i = 0; i < nx + 2; i++)
    {
        mat[i][0] = mat[i][1];
        mat[i][nx + 1] = mat[i][nx];
    }

    for (int i = 0; i < nx + 2; i++)
    {
        for (int j = 0; j < ny + 2; j++)
        {
            int ij = i + j * (nx + 2);
            charge_vec[ij] = mat[i][j];
        }
    }

    // update field Charge
    charge->val = charge_vec;
}

void species::advance_cell(int counter)
{
    float posx = vec.get()[counter].x;
    float posy = vec.get()[counter].y;

    float ux = vec.get()[counter].ux;
    float uy = vec.get()[counter].uy;
    float uz = vec.get()[counter].uz;

    int checker = 0;
    if (posx < 0.f)
    {
        vec.get()[counter].x = dx - posx;
        vec.get()[counter].ix = vec.get()[counter].ix - 1;
    }

    if (posx > dx)
    {
        vec.get()[counter].x = posx - dx;
        vec.get()[counter].iy = vec.get()[counter].ix + 1;
    }

    if (posy < 0.f)
    {
        vec.get()[counter].y = dy - posy;
        vec.get()[counter].iy = vec.get()[counter].iy - 1;
    }

    if (posy > dy)
    {
        vec.get()[counter].x = posy - dy;
        vec.get()[counter].iy = vec.get()[counter].iy + 1;
    }

    //! check boundaries to export data to mpi: later
}

void species::get_grid_points(int &nx_a, int &ny_a)
{
    nx_a = nx;
    ny_a = ny;
}

void species::print()
{
    for (int i = 0; i < np; i++)
    {
        std::cout << "cell (" << vec.get()[i].ix << ", " << vec.get()[i].iy << ")" << std::endl;
        std::cout << "x:" << vec.get()[i].x << std::endl;
        std::cout << "y: " << vec.get()[i].y << std::endl;
        std::cout << "ux: " << vec.get()[i].ux << std::endl;
        std::cout << "uy: " << vec.get()[i].uy << std::endl;
        std::cout << "****************" << std::endl;
    }
}

void species::write_output(int rank, int time)
{
    std::fstream Output_file;
    std::string filename;

    filename = "../results/grid_" + std::to_string(rank) + "__t_" + std::to_string(time) + ".txt";
    std::string space = "   ";

    Output_file.open(filename, std::ios::out);
    Output_file << "**************Grid**************" << std::endl;

    int precision = 4;
    for (int i = 0; i < np; i++)
    {
        Output_file << "cell (" << vec.get()[i].ix << ", " << vec.get()[i].iy << ")" << space;
        Output_file << "x:" << std::setw(precision) << vec.get()[i].x << space;
        Output_file << "y: " << std::setw(precision) << vec.get()[i].y << space;
        Output_file << "ux: " << std::setw(precision) << vec.get()[i].ux << space;
        Output_file << "uy: " << std::setw(precision) << vec.get()[i].uy << space;
        Output_file << std::endl;
    }
    Output_file << std::endl
                << std::endl;
    for (int i = 0; i < ny + 2; i++)
    {
        for (int j = 0; j < nx + 2; j++)
        {
            Output_file << std::setw(precision) << mat[i][j] << space;
        }
        Output_file << std::endl;
    }

    Output_file.close();
}