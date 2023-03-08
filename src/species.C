#include "species.h"

using namespace std;

// declare all stuff
species::species(string name_a, int *ppc_a, int *range_a, int *vec_U_a, float *vf_a, float *vth_a)
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
    np = set_nb_part();

    // initializing the array of particles
    vec = unique_ptr<part>(new part[np]);

    // random number generator
    random_device dev;
    mt19937_64 rng_aux(dev());
    normal_distribution<double> norm_aux(0, 1);
    rng = rng_aux;
    rand_gauss = norm_aux;
}

species::~species()
{
    cout << "destructor" << endl;
}

int species::set_nb_part()
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
    vector<float> loccell;

    const int npcell = ppc[0] * ppc[1];
    const float dpcellx = 1.0 / ppc[0];
    const float dpcelly = 1.0 / ppc[1];

    loccell.reserve(np);

    for (int j = 0; j < ppc[1]; j++)
    {
        for (int i = 0; i < ppc[0]; i++)
        {
            loccell.push_back(dpcellx * ((float)i + 0.5f)); // define axes inside grid cell
            loccell.push_back(dpcelly * ((float)j + 0.5f));
        }
    }

    for (int i = 0; i < loccell.size(); i++)
    {
        cout << "i: " << i << " value: " << loccell[i] << endl;
    }

    int ip = 0;
    // Uniform Density of Particles
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
}

void species::get_current(vector<float> &current_vec)
{
    current_vec.assign((range[1] * range[3]), 0.);

    int index = 0;
    for (int i = 0; i < np; i++)
    {
        index = vec.get()[i].ix + range[1] * vec.get()[i].iy;
        current_vec[index] = current_vec[index] + 1.;
    }
}

void species::advance()
{
}

void species::print()
{
    for (int i = 0; i < np; i++)
    {
        cout << "(" << vec.get()[i].ix << ", " << vec.get()[i].iy << ")" << endl;
        cout << "x:" << vec.get()[i].x << endl;
        cout << "y: " << vec.get()[i].y << endl;
        cout << "ux: " << vec.get()[i].ux << endl;
        cout << "uy: " << vec.get()[i].uy << endl;
        cout << "****************" << endl;
    }
}