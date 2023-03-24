#include "simulation.hh"
#include "species.hh"

void print(species *test, const int i)
{
    std::cout << "cell (" << test->vec[i].ix << ", " << test->vec[i].iy << ")" << std::endl;
    std::cout << "x:" << test->vec[i].x << std::endl;
    std::cout << "y: " << test->vec[i].y << std::endl;
    std::cout << "ux: " << test->vec[i].ux << std::endl;
    std::cout << "uy: " << test->vec[i].uy << std::endl;
    // std::cout << "++++++" << std::endl;
}

#define LOG(x)                       \
    {                                \
        std::cout << x << std::endl; \
    }

int main()
{
    std::string name = "electron";

    int ppc[2] = {1, 1};
    int range[2] = {10, 10}; // number of cells in each direction

    float *vf = new float[3];
    float vth[3] = {0.0, 0.0, 0.0};

    vf[0] = 0;
    vf[1] = 1.3;
    vf[2] = 0;

    float Ex = 0.;
    float Ey = 0.;
    int counter = 1;
    int flag_mpi = -2;

    species test(name, ppc, range, vf, vth);
    test.set_x();
    test.set_u();

    test.init_pusher(Ex, Ey, counter);

    for (int i = 0; i < 40; i++)
    {
        // std::cout << "pusher " << std::to_string(i) << std::endl;
        for (counter = 0; counter < test.np; counter++)
        {
            // std::cout << "pusher " << std::to_string(i) << std::endl;
            // print(&test, counter);
            test.particle_pusher(Ex, Ey, counter);
            // print(&test, counter);
            // std::cout << "advance" << std::endl;
            test.advance_cell(counter, flag_mpi);
            test.prepare_to_buffer(flag_mpi);
            // print(&test, counter);
            // std::cout << "****************" << std::endl;
        }
        test.write_output_vec(15, i);
    }

    // //Debugging buffer exchange information
    // LOG("here");
    // std::cout << "buffer south" << std::endl;
    // for (int i = 0; i < test.send_buffer_south.size(); i++)
    // {
    //     std::cout << "ix: " << test.send_buffer_south[i].ix << std::endl;
    //     std::cout << "iy: " << test.send_buffer_south[i].iy << std::endl;
    //     std::cout << "x: " << test.send_buffer_south[i].x << std::endl;
    //     std::cout << "y: " << test.send_buffer_south[i].y << std::endl;
    //     std::cout << "ux: " << test.send_buffer_south[i].ux << std::endl;
    //     std::cout << "uy: " << test.send_buffer_south[i].uy << std::endl;
    //     std::cout << "****************" << std::endl;
    // }

    // std::cout << "buffer north" << std::endl;
    // for (int i = 0; i < test.send_buffer_north.size(); i++)
    // {
    //     std::cout << "ix: " << test.send_buffer_north[i].ix << std::endl;
    //     std::cout << "iy: " << test.send_buffer_north[i].iy << std::endl;
    //     std::cout << "x: " << test.send_buffer_north[i].x << std::endl;
    //     std::cout << "y: " << test.send_buffer_north[i].y << std::endl;
    //     std::cout << "ux: " << test.send_buffer_north[i].ux << std::endl;
    //     std::cout << "uy: " << test.send_buffer_north[i].uy << std::endl;
    //     std::cout << "****************" << std::endl;
    // }

    // std::cout << "buffer east" << std::endl;
    // for (int i = 0; i < test.send_buffer_east.size(); i++)
    // {
    //     std::cout << "ix: " << test.send_buffer_east[i].ix << std::endl;
    //     std::cout << "iy: " << test.send_buffer_east[i].iy << std::endl;
    //     std::cout << "x: " << test.send_buffer_east[i].x << std::endl;
    //     std::cout << "y: " << test.send_buffer_east[i].y << std::endl;
    //     std::cout << "ux: " << test.send_buffer_east[i].ux << std::endl;
    //     std::cout << "uy: " << test.send_buffer_east[i].uy << std::endl;
    //     std::cout << "****************" << std::endl;
    // }

    // std::cout << "buffer west" << std::endl;
    // for (int i = 0; i < test.send_buffer_west.size(); i++)
    // {
    //     std::cout << "ix: " << test.send_buffer_west[i].ix << std::endl;
    //     std::cout << "iy: " << test.send_buffer_west[i].iy << std::endl;
    //     std::cout << "x: " << test.send_buffer_west[i].x << std::endl;
    //     std::cout << "y: " << test.send_buffer_west[i].y << std::endl;
    //     std::cout << "ux: " << test.send_buffer_west[i].ux << std::endl;
    //     std::cout << "uy: " << test.send_buffer_west[i].uy << std::endl;
    //     std::cout << "****************" << std::endl;
    // }

    // std::cout << "buffer ne" << std::endl;
    // for (int i = 0; i < test.send_buffer_ne.size(); i++)
    // {
    //     std::cout << "ix: " << test.send_buffer_ne[i].ix << std::endl;
    //     std::cout << "iy: " << test.send_buffer_ne[i].iy << std::endl;
    //     std::cout << "x: " << test.send_buffer_ne[i].x << std::endl;
    //     std::cout << "y: " << test.send_buffer_ne[i].y << std::endl;
    //     std::cout << "ux: " << test.send_buffer_ne[i].ux << std::endl;
    //     std::cout << "uy: " << test.send_buffer_ne[i].uy << std::endl;
    //     std::cout << "****************" << std::endl;
    // }

    // std::cout << "buffer nw" << std::endl;
    // for (int i = 0; i < test.send_buffer_nw.size(); i++)
    // {
    //     std::cout << "ix: " << test.send_buffer_nw[i].ix << std::endl;
    //     std::cout << "iy: " << test.send_buffer_nw[i].iy << std::endl;
    //     std::cout << "x: " << test.send_buffer_nw[i].x << std::endl;
    //     std::cout << "y: " << test.send_buffer_nw[i].y << std::endl;
    //     std::cout << "ux: " << test.send_buffer_nw[i].ux << std::endl;
    //     std::cout << "uy: " << test.send_buffer_nw[i].uy << std::endl;
    //     std::cout << "****************" << std::endl;
    // }

    // std::cout << "buffer sw" << std::endl;
    // for (int i = 0; i < test.send_buffer_sw.size(); i++)
    // {
    //     std::cout << "ix: " << test.send_buffer_sw[i].ix << std::endl;
    //     std::cout << "iy: " << test.send_buffer_sw[i].iy << std::endl;
    //     std::cout << "x: " << test.send_buffer_sw[i].x << std::endl;
    //     std::cout << "y: " << test.send_buffer_sw[i].y << std::endl;
    //     std::cout << "ux: " << test.send_buffer_sw[i].ux << std::endl;
    //     std::cout << "uy: " << test.send_buffer_sw[i].uy << std::endl;
    //     std::cout << "****************" << std::endl;
    // }

    // std::cout << "buffer se" << std::endl;
    // for (int i = 0; i < test.send_buffer_se.size(); i++)
    // {
    //     std::cout << "ix: " << test.send_buffer_se[i].ix << std::endl;
    //     std::cout << "iy: " << test.send_buffer_se[i].iy << std::endl;
    //     std::cout << "x: " << test.send_buffer_se[i].x << std::endl;
    //     std::cout << "y: " << test.send_buffer_se[i].y << std::endl;
    //     std::cout << "ux: " << test.send_buffer_se[i].ux << std::endl;
    //     std::cout << "uy: " << test.send_buffer_se[i].uy << std::endl;
    //     std::cout << "****************" << std::endl;
    // }

    test.write_output_buffer(20, 10);
    // LOG("160");
    test.get_charge();
    // LOG("charge");
    test.write_output_vec(20, 10);

    // test.update_part_list(); //with no sense here

    return 0;
}