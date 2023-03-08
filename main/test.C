#include "species.h"
#include <iostream>

using namespace std;

int main()
{
    string name = "electron";

    int ppc[2] = {1, 1};
    int range[4] = {0, 5, 0, 5}; // the latter ones are not included
    int vec_u[2] = {0, 20};

    float vf[3] = {0.1, 0.1, 0.1};
    float vth[3] = {0.1, 0.1, 0.1};

    species test(name, ppc, range, vec_u, vf, vth);
    test.set_X();
    test.set_U();
    cout << test.set_nb_part() << endl;
    cout << "********" << endl;
    // test.print();

    vector<float> charge;
    test.get_charge(charge);

    for (int i = 0; i < charge.size(); i++)
    {
        cout << "i: " << i << " I: " << charge[i] << endl;
    }

    return 0;
}
