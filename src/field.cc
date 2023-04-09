// #include "simulation.hh"
#include "field.hh"
#include <stdexcept>

namespace FCPIC
{
    // Constructors

    field::field(FCPIC_base const *base) : FCPIC_base(*base)
    {
        val.assign(N, 0.);
    }

    field::field(double *arr, FCPIC_base const *base) : FCPIC_base(*base)
    {
        val.assign(N, 0.);
        for (int i = 0; i < N; i++)
            val[i] = arr[i];
    }

    field::field(std::vector<double> &arr, FCPIC_base const *base) : FCPIC_base(*base)
    {
        val.reserve(N);
        val = arr;
    }

    void field::add_field(field *chargeA)
    {
        if ((N_x == chargeA->N_x) && (N_y == chargeA->N_y))
        {
            for (int i = 0; i < N; i++)
                val[i] = val[i] + chargeA->val[i];
        }
        else
        {
            throw std::invalid_argument("field - out of bound - check dimensions");
        }
    }

    field::~field()
    {
        val.clear();
    }

    void field::setValue(double value)
    {
        for (int i = 0; i < N; i++)
            val[i] = value;
    }

    void field::setNorthGuard(double *arr)
    {
        for (int j = 0; j < N_x; j++)
            val[NORTH_GUARD] = arr[j];
    }

    void field::setSouthGuard(double *arr)
    {
        for (int j = 0; j < N_x; j++)
            val[SOUTH_GUARD] = arr[j];
    }

    void field::setWestGuard(double *arr)
    {
        for (int i = 1; i <= N_int_y; i++)
            val[WEST_GUARD] = arr[i - 1];
    }

    void field::setEastGuard(double *arr)
    {
        for (int i = 1; i <= N_int_y; i++)
            val[EAST_GUARD] = arr[i - 1];
    }

    void field::reduceNorthBound(double *arr)
    {
        for (int j = 0; j < N_x; j++)
            val[NORTH_BOUND] += arr[j];
    }

    void field::reduceSouthBound(double *arr)
    {
        for (int j = 0; j < N_x; j++)
            val[SOUTH_BOUND] += arr[j];
    }

    void field::reduceWestBound(double *arr)
    {
        for (int i = 1; i <= N_int_y; i++)
            val[WEST_BOUND] += arr[i - 1];
    }

    void field::reduceEastBound(double *arr)
    {
        for (int i = 1; i <= N_int_y; i++)
            val[EAST_BOUND] += arr[i - 1];
    }

    void field::print_field(std::ostream &stream)
    {
        for (int i = N_y - 1; i >= 0; i--)
        {
            for (int j = 0; j < N_x; j++)
            {
                if (i == 0 || i == N_y - 1 || j == 0 || j == N_x - 1)
                    stream << " (" << this->val[POSITION] << ") ";
                else
                    stream << "  " << this->val[POSITION] << "  ";
            }
            stream << "\n";
        }
    }

    // change function
    // void field::field_grad()
    // {
    //     int nx_bulk = nx - 1;
    //     int ny_bulk = ny - 1;

    //     int idx = 0;
    //     float Ex = 0;
    //     float Ey = 0;

    //     // E calculation inside the domain
    //     for (int i = 1; i < nx_bulk; i++)
    //     {
    //         for (int j = 1; j < ny_bulk; j++)
    //         {
    //             idx = i + j * nx;
    //             Ex = -(pot_vec[idx + 1] - pot_vec[idx - 1]) / (2.f * dx);
    //             Ey = -(pot_vec[idx + nx] - pot_vec[idx - nx]) / (2.f * dy);
    //             Ex_vec[idx] = Ex;
    //             Ey_vec[idx] = Ey;
    //         }
    //     }
    //     // E calculation at the boundary conditions
    //     for (int i = 0; i < nx; i++)
    //     {
    //         int idx0 = i * nx;
    //         Ex_vec[idx0] = -(pot_vec[1 + idx0] - pot_vec[idx0]) / (2.f * dx);
    //         int idx1 = nx - 1 + i * ny;
    //         Ex_vec[idx1] = -(pot_vec[idx1] - pot_vec[idx1 - 1]) / (2.f * dx);
    //     }

    //     for (int i = 0; i < ny; i++)
    //     {
    //         int idy0 = i;
    //         Ey_vec[idy0] = -(pot_vec[nx + idy0] - pot_vec[idy0]) / (2.f * dy);
    //         int idy1 = i + (ny - 1) * ny;
    //         Ey_vec[idy1] = -(pot_vec[idy1] - pot_vec[idy1 - nx]) / (2.f * dy);
    //     }
    // }

    // void field::get_Efield(part &A, float &Ex_inter, float &Ey_inter)
    // {
    //     int icell = A.ix;
    //     int jcell = A.iy;
    //     float wx = A.x;
    //     float wy = A.y;

    //     int ij = A.ix + nx * A.iy;
    //     cout << "ij: " << ij << endl;
    //     cout << "x: " << wx << " y: " << wy << endl;

    //     float Aij = (dx - wx) * (dy - wy);
    //     float Aiij = wx * (dy - wy);
    //     float Aijj = (dx - wx) * wy;
    //     float Aiijj = wx * wy;

    //     cout << Aij << " " << Aiij << " " << Aijj << " " << Aiijj << endl;
    //     cout << Ex_vec[ij] << " " << Ex_vec[ij + 1] << " " << Ex_vec[ij + nx] << " " << Ex_vec[ij + nx + 1] << endl;
    //     cout << Ey_vec[ij] << " " << Ey_vec[ij + 1] << " " << Ey_vec[ij + nx] << " " << Ey_vec[ij + nx + 1] << endl;

    //     Ex_inter = Ex_vec[ij] * Aiij + Ex_vec[ij + 1] * Aiij + Ex_vec[ij + nx] * Aijj + Ex_vec[ij + 1 + nx] / (dx * dy);
    //     Ey_inter = Ey_vec[ij] * Aiij + Ey_vec[ij + 1] * Aiij + Ey_vec[ij + nx] * Aijj + Ey_vec[ij + 1 + nx] / (dx * dy);
    // }
}
