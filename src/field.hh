#ifndef _FIELD_
#define _FIELD_

#include <iostream>
#include <vector>

// Cycle in j from 0 to N_x (NOT INCLUDING) for(int j=0; j<N_x; j++)
#define SOUTH_GUARD j
#define NORTH_GUARD N_x *(N_y - 1) + j

// Cycle in i from 1 to N_int_y (INCLUDING) for(int i=1; i<=N_int_y; i++)
#define WEST_GUARD N_x *i
#define EAST_GUARD N_x *i + N_x - 1

// Cycle in j from 0 to N_x (NOT INCLUDING) for(int j=0; j<N_x; j++)
#define SOUTH_BOUND j + N_x
#define NORTH_BOUND N_x *(N_y - 2) + j

// Cycle in i from 1 to N_int_y (INCLUDING) for(int i=1; i<=N_int_y; i++)
#define WEST_BOUND N_x *i + 1
#define EAST_BOUND N_x *i + N_x - 2

// To range over internal cells:
// Cycle in i from 1 to N_int_y (INCLUDING) for(int i=1; i<=N_int_y; i++)
// Cycle in j from 1 to N_int_x (INCLUDING) for(int j=1; j<=N_int_x; j++)
// To range over all cells, including guard:
// Cycle in i from 0 to N_y (NOT INCLUDING) for(int i=0; i<N_y; i++)
// Cycle in j from 0 to N_x (NOT INCLUDING) for(int j=0; j<N_x; j++)
#define POSITION i *N_x + j
#define SOUTH i *N_x + j - N_x
#define NORTH i *N_x + j + N_x
#define WEST i *N_x + j - 1
#define EAST i *N_x + j + 1

namespace FCPIC
{
    // Boundary condition type
    typedef enum
    {
        TBD,
        PERIODIC,
        CONDUCTIVE,
    } BC_type;

    // Direction
    typedef enum
    {
        X_DIR,
        Y_DIR,
    } Direction;

    // Direction
    typedef enum
    {
        REPLACE,
        ADD,
    } Operation;

    // Field class
    class field
    {
    public:
        // Constructors
        // allocates memory to the field variables equal to the number of cells in the domain
        field();
        field(int, int);
        field(int, int, double *);
        field(int, int, std::vector<double> &);

        ~field(); // Destructor

        void setValue(double);

        void setNorthGuard(double *);
        void setSouthGuard(double *);
        void setEastGuard(double *);
        void setWestGuard(double *);

        void reduceNorthBound(double *);
        void reduceSouthBound(double *);
        void reduceWestBound(double *);
        void reduceEastBound(double *);

        void print_field(std::ostream &);

        // Member functions
        // void field_grad();                         // Uses central differences to calculate E field from Potential
        // void get_Efield(part &, float &, float &); // Interpolates the field inside the cell at particles' position

        std::vector<double> val;

        // worth putting public instead of get and set stuff
        int N;        // total size x*y
        int N_x, N_y; // total size

    private:
        // Members
        int N_int_x, N_int_y; // inner size
    };
}

#endif
