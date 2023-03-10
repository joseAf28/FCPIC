#include "fields.h"

using namespace std;

fields::fields(const int nx_a, const int ny_a, vector<float> &charge_a)
{
    nx = nx_a;
    ny = ny_a;

    charge_vec.reserve(nx * ny);

    int counter = nx;

    //! remove guard cell: it could be changed later
    //! not very efficient: copying the array again - change later
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
    int nx_bulk = nx - 1;
    int ny_bulk = ny - 1;

    int idx = 0;
    float Ex = 0;
    float Ey = 0;

    // E calculation inside the domain
    for (int i = 1; i < nx_bulk; i++)
    {
        for (int j = 1; j < ny_bulk; j++)
        {
            idx = i + j * nx;
            Ex = -(pot_vec[idx + 1] - pot_vec[idx - 1]) / (2.f * dx);
            Ey = -(pot_vec[idx + nx] - pot_vec[idx - nx]) / (2.f * dy);
            Ex_vec[idx] = Ex;
            Ey_vec[idx] = Ey;
        }
    }
    // E calculation at the boundary conditions
    for (int i = 0; i < nx; i++)
    {
        int idx0 = i * nx;
        Ex_vec[idx0] = -(pot_vec[1 + idx0] - pot_vec[idx0]) / (2.f * dx);
        int idx1 = nx - 1 + i * ny;
        Ex_vec[idx1] = -(pot_vec[idx1] - pot_vec[idx1 - 1]) / (2.f * dx);
    }

    for (int i = 0; i < ny; i++)
    {
        int idy0 = i;
        Ey_vec[idy0] = -(pot_vec[nx + idy0] - pot_vec[idy0]) / (2.f * dy);
        int idy1 = i + (ny - 1) * ny;
        Ey_vec[idy1] = -(pot_vec[idy1] - pot_vec[idy1 - nx]) / (2.f * dy);
    }
}

void fields::field_inter(part &A, float &Ex_inter, float &Ey_inter)
{
    int icell = A.ix;
    int jcell = A.iy;
    float wx = A.x;
    float wy = A.y;

    int ij = A.ix + nx * A.iy;
    cout << "ij: " << ij << endl;
    cout << "x: " << wx << " y: " << wy << endl;

    float Aij = (dx - wx) * (dy - wy);
    float Aiij = wx * (dy - wy);
    float Aijj = (dx - wx) * wy;
    float Aiijj = wx * wy;

    cout << Aij << " " << Aiij << " " << Aijj << " " << Aiijj << endl;
    cout << Ex_vec[ij] << " " << Ex_vec[ij + 1] << " " << Ex_vec[ij + nx] << " " << Ex_vec[ij + nx + 1] << endl;
    cout << Ey_vec[ij] << " " << Ey_vec[ij + 1] << " " << Ey_vec[ij + nx] << " " << Ey_vec[ij + nx + 1] << endl;

    Ex_inter = Ex_vec[ij] * Aiij + Ex_vec[ij + 1] * Aiij + Ex_vec[ij + nx] * Aijj + Ex_vec[ij + 1 + nx] / (dx * dy);
    Ey_inter = Ey_vec[ij] * Aiij + Ey_vec[ij + 1] * Aiij + Ey_vec[ij + nx] * Aijj + Ey_vec[ij + 1 + nx] / (dx * dy);
}

void fields::print()
{
    cout << "Potential Grid" << endl;
    cout << "**************" << endl;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            int idx = j + i * nx;
            cout << fixed << setprecision(4) << pot_vec[idx] << setw(5) << "  ";
        }
        cout << endl;
    }
    cout << "**************" << endl;

    cout << "Ex Grid" << endl;
    cout << "**************" << endl;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            int idx = j + i * nx;
            cout << fixed << setprecision(4) << Ex_vec[idx] << setw(5) << "  ";
        }
        cout << endl;
    }
    cout << "**************" << endl;
    cout << "Ey Grid" << endl;
    cout << "**************" << endl;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            int idx = j + i * nx;
            cout << fixed << setprecision(4) << Ey_vec[idx] << setw(5) << "  ";
        }
        cout << endl;
    }
    cout << "**************" << endl;
}