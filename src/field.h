#ifndef _FIELD_
#define _FIELD_

#include <iostream>
#include <vector>

namespace simulation
{
    // Boundary condition type
    typedef enum
    {
        NONE,
        PHYSICAL,
        BUFFER,
    } BC_type;

    // Direction
    typedef enum
    {
        X_DIR,
        Y_DIR,
    } Direction;

    // Field class
    class field
    {
    public:
        // Constructors
        // allocates memory to the field variables equal to the number of cells in the domain
        field(int, int);
        field(int, int, std::vector<double> &);
        field(int, int, float **);
        field();
        field(const field &);

        ~field(); // Destructor

        // Member functions
        void set_field_value(double);
        // void field_grad();                         // Uses central differences to calculate E field from Potential
        // void get_Efield(part &, float &, float &); // Interpolates the field inside the cell at particles' position

        // Members
        int Nx, Ny; // size
        int N;      // size
        std::vector<double> val;
        BC_type *bc;      // BC type
        double bc_val[5]; // BC value
        double h = 0.01;

        std::vector<float> Ex_vec;
        std::vector<float> Ey_vec;
    };
}

#endif
