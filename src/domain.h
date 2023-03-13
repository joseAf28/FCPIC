#ifndef _DOMAIN_
#define _DOMAIN_

#include "field.h"

namespace simulation
{
    class domain
    {
    public:
        // Members
        field *u;
        field *charge;

        // Constructor
        // adds the reference to the E and charge fields to the domain
        domain(field *Efield, field *con) : u(Efield), charge(con)
        {
            std::cout << __PRETTY_FUNCTION__ << std::endl;
        }
    };
}

#endif