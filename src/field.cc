// #include "simulation.hh"
#include "field.hh"
#include <stdexcept>

namespace FCPIC
{
    // Constructors
    // allocates memory to the field variables equal to the number of cells in the domain
    field::field()
    {
    }

    field::field(int Nx, int Ny)
    {
        N_int_x = Nx;
        N_int_y = Ny;
        N_x = Nx + 2;
        N_y = Ny + 2;
        N = N_x * N_y;
        val.assign(N, 0.);
        // std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    field::field(int Nx, int Ny, double *arr)
    {
        N_int_x = Nx;
        N_int_y = Ny;
        N_x = Nx + 2;
        N_y = Ny + 2;
        N = N_x * N_y;
        val.assign(N, 0.);
        for (int i = 0; i < N; i++)
            val[i] = arr[i];
    }

    field::field(int Nx, int Ny, std::vector<double> &arr)
    {
        N_int_x = Nx;
        N_int_y = Ny;
        N_x = Nx + 2;
        N_y = Ny + 2;
        N = N_x * N_y;
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
}
