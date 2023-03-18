#ifndef _FIELD_
#define _FIELD_

#include <iostream>
#include <vector>

//Cycle in j from 0 to N_x (NOT INCLUDING) for(int j=0; j<N_x; j++)
#define NORTH_GUARD j
#define SOUTH_GUARD (N_x-1)*N_y + j+1

//Cycle in i from 1 to N_int_y (INCLUDING) for(int i=1; i<=N_int_y; i++)
#define WEST_GUARD N_x*i 
#define EAST_GUARD 2*N_x*i -1

//Cycle in j from 0 to N_x (NOT INCLUDING) for(int j=0; j<N_x; j++)
#define NORTH_BOUND j+N_x
#define SOUTH_BOUND (N_x-2)*N_y + j+1

//Cycle in i from 1 to N_int_y (INCLUDING) for(int i=1; i<=N_int_y; i++)
#define WEST_BOUND N_x*i +1
#define EAST_BOUND 2*N_x*i -2

//Cycle in i from 1 to N_int_y (INCLUDING) for(int i=1; i<=N_int_y; i++)
//Cycle in j from 1 to N_int_x (INCLUDING) for(int j=1; j<=N_int_x; j++)
#define POSITION i*N_int_x + j
#define NORTH i*N_int_x + j - N_int_x
#define SOUTH i*N_int_x + j + N_int_x
#define WEST i*N_int_x + j - 1
#define EAST i*N_int_x + j + 1

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

    // Field class
    class field {
    public:
        // Constructors
        // allocates memory to the field variables equal to the number of cells in the domain
        field(int, int);
        //field(int, int, std::vector<double> &);
        //field(int, int, float **);
        //field();
        //field(const field &);

        ~field(); // Destructor

        void setValue(double);
    

        void getNorthBound(double*);
        void getSouthBound(double*);
        void getEastBound(double*);
        void getWestBound(double*);

        void setNorthGuard(double*);
        void setSouthGuard(double*);
        void setEastGuard(double*);
        void setWestGuard(double*);

        // Member functions
        // void field_grad();                         // Uses central differences to calculate E field from Potential
        // void get_Efield(part &, float &, float &); // Interpolates the field inside the cell at particles' position

        std::vector<double> val;
        
    private:
        // Members
        int N_x, N_y; // total size 
        int N;      // total size x*y
        int N_int_x, N_int_y; //inner size
    };
    }


#endif
