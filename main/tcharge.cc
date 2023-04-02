#include "simulation.hh"
#include "species.hh"
#include <unistd.h>

int main(int argc, char **argv)
{
    // not used for now
    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::normal_distribution<double> norm(3, 0.5);

    // initializaing simulation fields and MPI
    // FCPIC::simulation *sim = new FCPIC::simulation(argc, argv);
    // sim->set_conductive_field_bc();

    // declaring species object
    std::string name = "electron";
    int ppc[2] = {1, 1};
    int range[2] = {10, 20}; // number of cells in each direction

    float *vfa = new float[3];
    float *vfb = new float[3];
    float vth[3] = {2., 2., 2.};

    vfa[0] = 0.5;
    vfa[1] = 0.3;
    vfa[2] = 0;
    vfb[0] = 0.4;
    vfb[1] = 0.4;
    vfb[2] = 0;

    FCPIC::field *Ex = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *Ey = new FCPIC::field(range[0] + 1, range[1] + 1);
    FCPIC::field *charge;
    FCPIC::field *phi = new FCPIC::field(range[0] + 1, range[1] + 1);

    species specA(name, ppc, range, vfa, vth, 1.);
    species specB(name, ppc, range, vfa, vth, 1.);
    specA.set_x();
    specA.set_u();

    specB.set_x();
    specB.set_u();

    charge = new FCPIC::field(range[0] + 1, range[1] + 1);
    specA.get_charge(charge);
    charge->print_field(std::cout);
    std::cout << std::endl;
    specB.get_charge(charge);

    charge->print_field(std::cout);

    charge->setValue(0.f);
    charge->print_field(std::cout);

    return 0;
}