#include "fields.h"

using namespace std;

fields::fields(const int nx_a, const int ny_a, vector<float> &charge_a)
{
    nx = nx_a;
    ny = ny_a;

    charge_vec.reserve(nx * ny);

    int counter = nx;

    //! remove guard cell: it could be changed later
    for (int i = 0; i < charge_a.size(); i++)
    {
        if (i == counter)
        {
            counter = counter + nx + 1.;
            // cout << "counter " << counter << "  i: " << i << endl;
            continue;
        }

        charge_vec.push_back(charge_a[i]);
    }

    // debugging charge_vec
    // cout << "*************" << endl;
    // for (int j = 0; j < nx; j++)
    // {
    //     for (int i = 0; i < ny; i++)
    //     {
    //         int idx = i + j * nx;
    //         cout << charge_vec[idx] << " ";
    //     }
    //     cout << endl;
    // }

    // initializing fields: set to 0. as intial conditions;
    pot_vec.assign(nx * ny, 0.f);
    Ex_vec.assign(nx * ny, 0.f);
    Ey_vec.assign(nx * ny, 0.f);
}

fields::~fields()
{
    cout << __PRETTY_FUNCTION__ << endl;
    charge_vec.clear();
    pot_vec.clear();
    Ex_vec.clear();
    Ey_vec.clear();
}

void fields::potential_solver()
{
    int nx_bulk = nx - 1;
    int ny_bulk = ny - 1;

    float error_norm = 1.f;
    float pot_aux = 0.f;

    int iter = 0;

    //! we impose dirichelet conditions pot(0, j) = pot(nx-1, j) = 0
    //!                                 pot(i, 0) = pot(i, ny-1) = 0
    //! stopping criteria |xk+1 - xk | < tol

    while (error_norm > tol && iter < max_iter)
    {
        // cout << "iter: " << iter << endl;
        // cout << "error_norm: " << setprecision(10) << error_norm << endl;

        error_norm = 0.f;
        // iteration m
        for (int i = 1; i < nx_bulk; i++)
        {
            for (int j = 1; j < ny_bulk; j++)
            {
                int idx = i + j * nx;
                pot_aux = pot_vec[idx];
                pot_vec[idx] = (pot_vec[idx + 1] + pot_vec[idx - 1] + pot_vec[idx - nx] + pot_vec[idx + nx] - charge_vec[idx] * dx * dy) / 4.f;
                error_norm = error_norm + abs(pot_vec[idx] - pot_aux);
            }
        }
        iter = iter + 1;
    }
}

void fields::field_solver()
{
}

void fields::print()
{
    cout << "**************" << endl;
    cout << "Potential Grid" << endl;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            int idx = i + j * nx;
            cout << fixed << setprecision(4) << pot_vec[idx] << setw(5) << "  ";
        }
        cout << endl;
    }
}