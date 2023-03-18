#include "species.h"

int main()
{
    std::string name = "electron";

    int ppc[2] = {3, 3};
    int range[4] = {0, 100, 0, 100}; // the latter ones are not included
    int vec_u[2] = {0, 20};

    float vf[3] = {0.1, 0.1, 0.1};
    float vth[3] = {0.1, 0.1, 0.1};

    species test(name, ppc, range, vec_u, vf, vth);
    test.set_X();
    test.set_U();
    test.get_charge();

    test.write_output(1);
    // for (int i = 0; i < charge.size(); i++)
    // {
    //     cout << "i: " << i << " I: " << charge[i] << endl;
    // }

    // box grid size info
    int nx = 0;
    int ny = 0;

    test.get_grid_points(nx, ny);
    std::cout << "nx: " << nx << "  ny: " << ny << std::endl;
    return 0;
}
