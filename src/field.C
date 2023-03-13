#include "simulation.h"
#include "field.h"

namespace simulation
{
    // Constructors
    // allocates memory to the field variables equal to the number of cells in the domain
    field::field(int N_x, int N_y)
    {
        Nx = N_x;
        Ny = N_y;
        N = N_x * N_y;
        val.assign(N, 0.);
        bc = new BC_type[N];
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    field::field(int N_x, int N_y, std::vector<double> &value)
    {
        Nx = N_x;
        Ny = N_y;
        N = N_x * N_y;
        val.reserve(N);
        val = value;
        bc = new BC_type[N];
    }

    field::field()
    {
        Nx = 0;
        Ny = 0;
        N = Nx * Ny;
        val.assign(N, 0.);
        bc = new BC_type[N];
    }

    field::field(const field &obj)
    {
        Nx = obj.Nx;
        Ny = obj.Ny;
        N = obj.N;
        val = obj.val;
        bc = new BC_type[N];
        memcpy(bc, obj.bc, sizeof(BC_type) * N);
    }

    void field::set_field_value(double value)
    {
        for (int i = 0; i < N; i++)
            this->val[i] = value;
    }

    field::~field()
    {
        val.clear();
        delete bc;
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
