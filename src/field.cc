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

    field::field(float *arr, FCPIC_base const *base) : FCPIC_base(*base)
    {
        val.assign(N, 0.);
        for (int i = 0; i < N; i++)
            val[i] = arr[i];
    }

    field::field(std::vector<float> &arr, FCPIC_base const *base) : FCPIC_base(*base)
    {
        val.reserve(N);
        val = arr;
    }

    field::~field()
    {
        val.clear();
    }

    void field::setValue(float value)
    {
        for (int i = 0; i < N; i++)
            val[i] = value;
    }

    void field::setNorthGuard(float *arr)
    {
        for (int j = 0; j < N_x; j++)
            val[NORTH_GUARD] = arr[j];
    }

    void field::setSouthGuard(float *arr)
    {
        for (int j = 0; j < N_x; j++)
            val[SOUTH_GUARD] = arr[j];
    }

    void field::setWestGuard(float *arr)
    {
        for (int i = 1; i <= N_int_y; i++)
            val[WEST_GUARD] = arr[i - 1];
    }

    void field::setEastGuard(float *arr)
    {
        for (int i = 1; i <= N_int_y; i++)
            val[EAST_GUARD] = arr[i - 1];
    }

    void field::reduceNorthBound(float *arr)
    {
        for (int j = 0; j < N_x; j++)
            val[NORTH_BOUND] += arr[j];
    }

    void field::reduceSouthBound(float *arr)
    {
        for (int j = 0; j < N_x; j++)
            val[SOUTH_BOUND] += arr[j];
    }

    void field::reduceWestBound(float *arr)
    {
        for (int i = 1; i <= N_int_y; i++)
            val[WEST_BOUND] += arr[i - 1];
    }

    void field::reduceEastBound(float *arr)
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
