//FCPIC - 2D Particle-in-Cell code using MPI
//Guilherme Crispim, João Palma, José Afonso
//Advanced Topics in Computational Physics, 2023, IST

//File field.cc:
//Implementation of all methods of class Field

#include "field.hh"
#include <stdexcept>

namespace FCPIC
{
    /* 
        Constructor of class Fields
        Inputs: FCPIC_base object
        + Creates a field grid with the dimensions specified in
        the FCPIC_base object, stored in a vector, with all 
        entries set to 0
    */
    field::field(FCPIC_base const *base) : FCPIC_base(*base)
    {
        val.assign(N, 0.);
    }

    /* 
        Constructor of class Fields
        Inputs: FCPIC_base object and float pointer
        + Creates a field grid with the dimensions specified in
        the FCPIC_base object, stored in a vector, with all entries 
        set to the values provided in the float array
    */
    field::field(float *arr, FCPIC_base const *base) : FCPIC_base(*base)
    {
        val.assign(N, 0.);
        for (int i = 0; i < N; i++)
            val[i] = arr[i];
    }

    /* 
        Constructor of class Fields
        Inputs: FCPIC_base object and vector of floats
        + Creates a field grid with the dimensions specified in
        the FCPIC_base object, stored in a vector, with all entries
        set to the values provided in the float vector 
    */
    field::field(std::vector<float> &arr, FCPIC_base const *base) : FCPIC_base(*base)
    {
        val.reserve(N);
        val = arr;
    }

    /* 
        Destructor of class Fields
        Clears the vector of field values
    */
    field::~field()
    {
        val.clear();
    }

    /* 
        Function setValue
        Inputs: float
        + Sets all field entries to the value provided
    */
    void field::setValue(float value)
    {
        for (int i = 0; i < N; i++)
            val[i] = value;
    }

    /* 
        Functions set[Direction]Guard
        Inputs: float pointer
        + Sets the guard cell values of the specified direction
        with the values provided in the float pointer
    */
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

    /* 
        Functions reduce[Direction]Bound
        Inputs: float pointer
        + Sums the values of the provided float pointer in
        the boundary cell values of the specified direction
    */
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

    /* 
        Function print_field
        Inputs: output stream for writing
        + Auxiliary function for printing the values of
        the grid
    */
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
