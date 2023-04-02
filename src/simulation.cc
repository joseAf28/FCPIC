#include "simulation.hh"
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <unistd.h>

#define epsilon0 55263494.06 //in e^2*eV^-1*m^-1

namespace FCPIC
{
    /********************************************************
    Solving the poisson equation in 2D parallely
    N_int_x and N_int_y are the number of inner grid points in x and y direction respectively
    ***********************************************************/

    simulation::simulation(int argc, char **argv)
    {
        int rank;
        MPI_Init(&argc, &argv);
        // retrieve the number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &n_Procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        readArgs(argc,argv);
        setParams();
        if(rank ==0)
            printTitle();

        
        aspect = 1; // (INPUT) y_len = aspect (x_len always norm to 1)
        N_int_x = 21; //
        N_int_y = 11; //
        N_x = N_int_x + 2;//
        N_y = N_int_y + 2;//
        N = N_x * N_y;//

        dx = 1; //
        dy = 1; //

        grid[X_DIR] = 2;//
        grid[Y_DIR] = 2;//

        Y_guard_data = new double[N_x];
        X_guard_data = new double[N_int_y];

        Y_guard_data1 = new double[N_x]();
        X_guard_data1 = new double[N_int_y]();
        Y_guard_data2 = new double[N_x]();
        X_guard_data2 = new double[N_int_y]();
    }

    simulation::~simulation()
    {
        MPI_Type_free(&exchange_field_type[X_DIR]);
        MPI_Type_free(&exchange_field_type[Y_DIR]);
        MPI_Type_free(&exchange_part_type);

        MPI_Finalize();

        delete[] Y_guard_data, X_guard_data,
            Y_guard_data1, X_guard_data1,
            Y_guard_data2, X_guard_data2;
    }

    std::string simulation::print_SI(double x){
        int exponent = floor(log10(x));
        exponent -= exponent%3;
        double mantissa = x/pow(10.,exponent);

        std::string output = std::to_string(mantissa);

        if(exponent == -12) output.append(" p");
        if(exponent == -9) output.append(" n");
        if(exponent == -6) output.append(" μ");
        if(exponent == -3) output.append(" m");
        if(exponent == 0) output.append(" ");
        if(exponent == 3) output.append(" k");
        if(exponent == 6) output.append(" M");
        if(exponent == 9) output.append(" G");
        if(exponent == 12) output.append(" T");
    
    return output;

    }

    void simulation::printTitle(){
        std::cout << "\n";
        std::cout << "┌─────────────────────────────────────────────────────────────────────────────┐\n";
        std::cout << "│ ▄▄▄▄▄▄▄▄▄▄▄▄   ▄▄▄▄▄▄▄▄▄▄▄▄   ▄▄▄▄▄▄▄▄▄▄▄▄▄   ▄▄▄▄▄▄▄▄▄▄▄▄▄   ▄▄▄▄▄▄▄▄▄▄▄▄▄ │\n";
        std::cout << "│ ████████████   ████████████   █████████████   █████████████   █████████████ │\n";
        std::cout << "│ ████████████   ████████████   █████████████   █████████████   █████████████ │\n";
        std::cout << "│ █████          █████          █████     ███       █████       █████         │\n";
        std::cout << "│ █████          █████          █████     ███       █████       █████         │\n";
        std::cout << "│ █████▄▄▄▄▄     █████          █████▄▄▄▄▄███       █████       █████         │\n";
        std::cout << "│ ██████████     █████          █████████████       █████       █████         │\n";
        std::cout << "│ █████          █████          █████               █████       █████         │\n";
        std::cout << "│ █████          █████▄▄▄▄▄▄▄   █████           ▄▄▄▄█████▄▄▄▄   █████▄▄▄▄▄▄▄▄ │\n";
        std::cout << "│ █████          ████████████   █████           █████████████   █████████████ │\n";
        std::cout << "│ █████          ████████████   █████           █████████████   █████████████ │\n";
        std::cout << "├─────────────────────────────────────────────────────────────────────────────┤\n";
        std::cout << "│ 2D Particle-in-Cell code using MPI                                          │\n";
        std::cout << "│ Guilherme Crispim, João Palma, José Afonso, ATCP 2023                       │\n";
        std::cout << "└─────────────────────────────────────────────────────────────────────────────┘\n";
        std::cout << "\n";
        std::cout << "Simulation size: " << print_SI(xlen*Lref) << "m x " << print_SI(xlen*aspect*Lref) << "m\n";
        std::cout << "Spatial discretization: "<< N_total_x << "x" << N_total_y << "\n";
        std::cout << "MPI process grid: " << n_Procs << " processes (" 
                  << grid[X_DIR]<< "x" << grid[Y_DIR] <<  ")\n";
        std::cout << "Electron Debye length: " << print_SI(Lref) << "m\n";
        std::cout << "Plasma frequency: " << print_SI(1/Tref) << "Hz\n";
        std::cout << "\n";
    }

    void simulation::printHelp(){
        std::cout << "\n";
        std::cout << ">> mpiexec ./FCPIC.exec -infile=infile.txt\n";
        std::cout << ">> mpiexec ./FCPIC.exec -npart=1000,2000 -charge=1,-1  ...\n";
        std::cout << "\n";
        std::cout << "-infile=fff           Read parameters from input file\n";
        std::cout << "\n";
        std::cout << "-npart=nn1,nn2,...    Particle number for all species\n";
        std::cout << "                            Default: 1000\n";
        std::cout << "-charge=qq1,qq2,...   Charge of all species (in q_proton/m_electron)\n";
        std::cout << "                            Default: -1\n";
        std::cout << "-mass=mm1,mm2,...     Mass of all species (in m_electron)\n";
        std::cout << "                            Default: 1\n";
        std::cout << "-temp=tt1,tt2,...     Temperature of all species (in eV)\n";
        std::cout << "                            Default: 1\n";
        std::cout << "-vxfluid=vx1,vx2,...  X component of fluid velocity (in v_thermal)\n";
        std::cout << "                            Default: 0\n";
        std::cout << "-vyfluid=vy1,vy2,...  Y component of fluid velocity (in v_thermal)\n";
        std::cout << "                            Default: 0\n";
        std::cout << "-nxlen=lll            Horizontal length of the simulation box (in m)\n";
        std::cout << "                            Default: 1\n";
        std::cout << "-nxproc=ppp           Number of MPI processes horizontally in the grid\n";
        std::cout << "                            Default: Number of processes/2\n";
        std::cout << "-nxsamples=sss        Number of discretization cells horizontally\n";
        std::cout << "                            Default: 100\n";
        std::cout << "-aspect=aaa           Box aspect ratio (ylen = aspect * xlen)\n";
        std::cout << "                            Default: 1\n";
        std::cout << "-simtime=ttt          Simulation time (in secs)\n";
        std::cout << "                            Default: 1\n";
        std::cout << "-boundcond=bbb        Boundary condition (1->periodic, 2->conductive)\n";
        std::cout << "                            Default: 0\n";
        std::cout << "\n"; 
    }

    void simulation::readArgs(int argc, char **argv){
        std::vector<std::string> allArgs(argv, argv+argc);
        std::vector<std::string> numbers;
        std::string line, header, number;
        std::vector<bool> def_values(12, true);
        int k;

        for(auto & arg : allArgs){
            k=0;
            for(char &c : arg){
                    if((c > 47 && c<58) || c == '.' || (c == '-' && k != 0))
                        number.push_back(c);
                    if(c == ','){
                        numbers.push_back(number);
                        number.clear();
                        }
                    if((c > 96 && c < 123) || (c == '-' && k == 0))
                        header.push_back(c);
                    if(c > 64 && c < 91)
                        header.push_back(c+32);

                    if(header.compare("-infile")==0){
                        arg.erase(0,8);
                        break;
                    }
                    k++;
                }

                if(header.compare("-help")==0){
                    printHelp();
                }

                if(header.compare("-infile")==0){
                    def_values.assign(12,true);
                    getParamsfromFile(arg, &def_values);
                    break;
                }

                if(header.compare("-npart")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        Npart.push_back(std::stoi(num));
                    def_values[0]=false;
                }
                if(header.compare("-charge")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        charge.push_back(std::stod(num));
                    def_values[1]=false;
                }
                if(header.compare("-mass")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        mass.push_back(std::stod(num));
                    def_values[2]=false;
                }
                if(header.compare("-temp")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        temp.push_back(std::stod(num));
                    def_values[3]=false;
                }
                if(header.compare("-vxfluid")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        vxfluid.push_back(std::stod(num));
                    def_values[4]=false;
                }
                if(header.compare("-vyfluid")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        vyfluid.push_back(std::stod(num));
                    def_values[5]=false;
                }
                if(header.compare("-xlen")==0){
                    xlen = stod(number);
                    def_values[6]=false;
                }
                if(header.compare("-nxproc")==0){
                    nxproc = stoi(number);
                    def_values[7]=false;
                }
                if(header.compare("-nxsamples")==0){
                    N_total_x = stoi(number);
                    def_values[8]=false;
                }
                if(header.compare("-aspect")==0){
                    aspect = stod(number);
                    def_values[9]=false;
                }
                if(header.compare("-simtime")==0){
                    simtime = stoi(number);
                    def_values[10]=false;
                }
                if(header.compare("-boundcond")==0){
                    bc[X_DIR] = stoi(number);
                    bc[Y_DIR] = stoi(number);
                    def_values[11]=false;
                }
                
                numbers.clear();
                number.clear();
                header.clear();
        }

        if(def_values[0]) Npart.push_back(1000);
        if(def_values[1]) charge.push_back(-1.);
        if(def_values[2]) mass.push_back(1.);
        if(def_values[3]) temp.push_back(1.);
        if(def_values[4]) vxfluid.push_back(0.);
        if(def_values[5]) vyfluid.push_back(0.);
        if(def_values[6]) xlen=1;
        if(def_values[7]) nxproc = n_Procs/2;
        if(def_values[8]) N_total_x = 100;
        if(def_values[9]) aspect = 1;
        if(def_values[10]) simtime = 1;
        if(def_values[11]){
            bc[X_DIR] = PERIODIC;
            bc[Y_DIR] = PERIODIC;
        }
    }

    void simulation::getParamsfromFile(std::string filename, std::vector<bool> *def_values){
        
            std::ifstream infile(filename);
            std::vector<std::string> filelines, numbers;
            std::string line, header, number;

            while(getline(infile, line))
                filelines.push_back(line);
            
            infile.close();
            
            for(auto & fline : filelines){
                for(char &c : fline){
                    if((c > 47 && c<58) || c == '.' || c == '-')
                        number.push_back(c);
                    if(c == ','){
                        numbers.push_back(number);
                        number.clear();
                        }
                    if(c > 96 && c < 123)
                        header.push_back(c);
                    if(c > 64 && c < 91)
                        header.push_back(c+32);
                    if(c == '#')
                        break;
                }
                
                if(header.compare("npart")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        Npart.push_back(std::stoi(num));
                    (*def_values)[0]=false;
                }
                if(header.compare("charge")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        charge.push_back(std::stod(num));
                    (*def_values)[1]=false;
                }
                if(header.compare("mass")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        mass.push_back(std::stod(num));
                    (*def_values)[2]=false;
                }
                if(header.compare("temp")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        temp.push_back(std::stod(num));
                    (*def_values)[3]=false;
                }
                if(header.compare("vxfluid")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        vxfluid.push_back(std::stod(num));
                    (*def_values)[4]=false;
                }
                if(header.compare("vyfluid")==0){
                    numbers.push_back(number);
                    for(auto & num : numbers)
                        vyfluid.push_back(std::stod(num));
                    (*def_values)[5]=false;
                }
                if(header.compare("xlen")==0){
                    xlen = stod(number);
                    (*def_values)[6]=false;
                }
                if(header.compare("nxproc")==0){
                    nxproc = stoi(number);
                    (*def_values)[7]=false;
                }
                if(header.compare("nxsamples")==0){
                    N_total_x = stoi(number)+1;
                    (*def_values)[8]=false;
                }
                if(header.compare("aspect")==0){
                    aspect = stod(number);
                    (*def_values)[9]=false;
                }
                if(header.compare("simtime")==0){
                    simtime = stoi(number);
                    (*def_values)[10]=false;
                }
                if(header.compare("boundcond")==0){
                    bc[X_DIR] = stoi(number);
                    bc[Y_DIR] = stoi(number);
                    (*def_values)[11]=false;
                }
                
                numbers.clear();
                number.clear();
                header.clear();
            }
    }

    void simulation::setParams(){
        if(Npart.size() != charge.size() || Npart.size() != mass.size() || 
           Npart.size() != temp.size() || Npart.size() != vxfluid.size() || 
           Npart.size() != vyfluid.size()){
            std::cout << "ERROR: PLEASE PROVIDE ARGS FOR ALL SPECIES\n";
        }

        Nspecies = Npart.size();

        for(int k = 0; k < Nspecies; k++){
            if(Npart[k] < 1){
                std::cout << "ERROR: NON-PHYSICAL NUMBER OF PARTICLES\n";
            }
            if(mass[k] <= 0){
                std::cout << "ERROR: NON-PHYSICAL TEMPERATURE\n";
            }
            if(temp[k] < 0){
                std::cout << "ERROR: NON-PHYSICAL TEMPERATURE\n";
            }
        }
        if(xlen <= 0){
            std::cout << "ERROR: NON-PHYSICAL BOX LENGTH\n";
        }
        if(n_Procs%nxproc != 0){
            std::cout << "ERROR: UNABLE TO FORM A PROCESS GRID\n";
        }
        if(N_total_x < 10){
            std::cout << "ERROR: NOT ENOUGH GRID CELLS\n";
        }
        if(aspect <= 0){
            std::cout << "ERROR: INVALID ASPECT RATIO\n";
        }
        if(simtime <= 0){
            std::cout << "ERROR: INVALID SIMULATION TIME\n";
        }
        if(bc[X_DIR] != 1 && bc[X_DIR] != 2){
            std::cout << "ERROR: INVALID BOUNDARY CONDITIONS\n";
        }

        //All norms are done in reference to the first species temperature and density
        Vref = 419382.88*sqrt(temp[0]); //Normalizing speeds to Vthermal of 1st species (m/s)

        for(int k = 0; k < Nspecies; k++){
                temp[k] /= Vref;
        }

        Nref = (double)Npart[0]/(aspect*xlen*xlen); //2D Density (m^-2)
        Lref = sqrt(55263494.06*temp[0]/pow(Nref,1.5)); //Normalizing lengths to electron Debye length (m)
        Tref = Lref/Vref; //Normalizing times to electron inverse plasma frequency (s)

        simtime /= Tref;
        xlen /= Lref;

        dx=1./3.; //spatial discretization must be ~ Debye length/3 for stability

        N_total_x = round((double)xlen/ (double)dx);
        N_total_y = round(aspect*(double)N_total_x);

        grid[X_DIR] = nxproc;
        grid[Y_DIR] = n_Procs/nxproc;

        if(bc[X_DIR] != PERIODIC)
            N_total_x -= 2;
        if(bc[X_DIR] != PERIODIC)
            N_total_y -= 2;

        int k = 0;
        while(N_total_x%grid[X_DIR] != 0){
            k = -k - abs(k)/k;
            N_total_x += k;
        }

        N_int_x = N_total_x/grid[X_DIR];            
        N_x = N_int_x + 2;

        k = 0;
        while(N_total_y%grid[Y_DIR] != 0){
            k = -k - abs(k)/k;
            N_total_y += k;
        }

        N_int_y = N_total_y/grid[Y_DIR];            
        N_y = N_int_y + 2;

        if(bc[X_DIR] != PERIODIC)
            N_total_x += 2;
        if(bc[X_DIR] != PERIODIC)
            N_total_y += 2;
        
        N = N_x * N_y;

        dx = xlen / (double)N_total_x;//
        dy = aspect*xlen / (double)N_total_y;//

        dt = 1/(std::max(1.,sqrt(vxfluid[0]*vxfluid[0]+vyfluid[0]*vyfluid[0]))*(1/dx + 1/dy));

        wrap_around[X_DIR] = bc[X_DIR]==PERIODIC ? 1 : 0;
        wrap_around[Y_DIR] = bc[Y_DIR]==PERIODIC ? 1 : 0;
    }

    // Creating a virtual cartesian topology
    void simulation::setup_proc_grid()
    {
        // number of processes per row and column

        int reorder = 1; // reorder process ranks

        // creates a new communicator, grid_comm
        MPI_Cart_create(MPI_COMM_WORLD, 2, grid, wrap_around, reorder, &grid_comm);

        // retrieve new rank and cartesian coordinates of this process
        MPI_Comm_rank(grid_comm, &grid_rank);
        MPI_Cart_coords(grid_comm, grid_rank, 2, grid_coord);

        // calculate ranks of neighboring processes in the grid
        MPI_Cart_shift(grid_comm, 1, 1, &grid_left, &grid_right);
        MPI_Cart_shift(grid_comm, 0, 1, &grid_bottom, &grid_top);

        // Datatype for Field's Communication
        // Datatype for horizontal data exchange
        MPI_Type_vector(N_int_y, 1, N_x, MPI_DOUBLE, &exchange_field_type[X_DIR]);
        MPI_Type_commit(&exchange_field_type[X_DIR]);

        // Datatype for vertical data exchange
        MPI_Type_vector(N_x, 1, 1, MPI_DOUBLE, &exchange_field_type[Y_DIR]);
        MPI_Type_commit(&exchange_field_type[Y_DIR]);

        // Datatype for Species's communication
        int blocklengths[8] = {1, 1, 1, 1, 1, 1, 1, 1};
        MPI_Datatype types[8] = {MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_CXX_BOOL};

        offsets[0] = offsetof(part, ix);
        offsets[1] = offsetof(part, iy);
        offsets[2] = offsetof(part, x);
        offsets[3] = offsetof(part, y);
        offsets[4] = offsetof(part, ux);
        offsets[5] = offsetof(part, uy);
        offsets[6] = offsetof(part, uz);
        offsets[7] = offsetof(part, flag);

        MPI_Type_create_struct(nitems, blocklengths, offsets, types, &exchange_part_type);
        MPI_Type_commit(&exchange_part_type);

        // obtaing the coordinates of the diagonal processes
        int coords_ne[2] = {grid_coord[X_DIR] + 1, grid_coord[Y_DIR] + 1};
        int coords_se[2] = {grid_coord[X_DIR] - 1, grid_coord[Y_DIR] + 1};
        int coords_nw[2] = {grid_coord[X_DIR] + 1, grid_coord[Y_DIR] - 1};
        int coords_sw[2] = {grid_coord[X_DIR] - 1, grid_coord[Y_DIR] - 1};

        get_diagonal_rank(coords_ne, grid_ne);
        get_diagonal_rank(coords_se, grid_se);
        get_diagonal_rank(coords_nw, grid_nw);
        get_diagonal_rank(coords_sw, grid_sw);

        //std::cout << "grid_rank: " << grid_rank << " ne: " << grid_ne << " se: " << grid_se << " sw: " << grid_sw << " nw: " << grid_nw << std::endl;

        /////////////////////777////////////

        //---TESTES---
        // PARA TESTAR COISAS USAR ESTE SÍTIO
        /*
        std::cout << "Grid rank: " << grid_rank << std::endl;
        std::cout << "N_int_x: " << N_int_x << "   N_int_y: " << N_int_y << std::endl;
        std::cout << "N_x: " << N_x << "   N_y: " << N_y << std::endl;
        double array[25] = { 0, 1, 2, 3, 4,
                             5, 6, 7, 8, 9,
                            10,11,12,13,14,
                            15,16,17,18,19,
                            20,21,22,23,24};
        for(int k = 0; k<25; k++)
            array[k] += grid_rank*25;
        field phi(N_int_x, N_int_y, array);
        //phi.setValue((double) grid_rank);
        sleep(grid_rank);
        phi.print_field(std::cout);
        exchange_charge_buffers(&phi);
        sleep(grid_rank);
        phi.print_field(std::cout);
        */
    }

    void simulation::get_diagonal_rank(int *coords, int &id_proc)
    {
        // both wrap around in X_DIR and Y_DIR must be 0 or 1 at the same time
        if (!wrap_around[X_DIR] && !wrap_around[Y_DIR] && (coords[0] >= grid[0] || coords[1] >= grid[1] || coords[0] < 0 || coords[1] < 0))
            id_proc = MPI_PROC_NULL;
        else
            MPI_Cart_rank(grid_comm, coords, &id_proc);
    }

    void simulation::set_periodic_field_bc()
    {
        bc[X_DIR] = PERIODIC;
        bc[Y_DIR] = PERIODIC;
        wrap_around[X_DIR] = 1;
        wrap_around[Y_DIR] = 1;

        setup_proc_grid();
    }

    void simulation::set_conductive_field_bc()
    {
        bc[X_DIR] = CONDUCTIVE;
        bc[Y_DIR] = CONDUCTIVE;
        wrap_around[X_DIR] = 0;
        wrap_around[Y_DIR] = 0;

        setup_proc_grid();
    }

    void simulation::exchange_phi_buffers(field *phi)
    {
        int i, j;
        
        i = 1;
        j = 0;

        // Communication stream leftward
        MPI_Sendrecv(&phi->val[WEST_BOUND], 1, exchange_field_type[X_DIR], grid_left, 0,
                     &phi->val[EAST_GUARD], 1, exchange_field_type[X_DIR], grid_right, 0,
                     grid_comm, &status);

        // Communication stream rightward
        MPI_Sendrecv(&phi->val[EAST_BOUND], 1, exchange_field_type[X_DIR], grid_right, 0,
                     &phi->val[WEST_GUARD], 1, exchange_field_type[X_DIR], grid_left, 0,
                     grid_comm, &status);

        // Communication stream upward
        MPI_Sendrecv(&phi->val[NORTH_BOUND], 1, exchange_field_type[Y_DIR], grid_top, 0,
                     &phi->val[SOUTH_GUARD], 1, exchange_field_type[Y_DIR], grid_bottom, 0,
                     grid_comm, &status);

        // Communication stream downward
        MPI_Sendrecv(&phi->val[SOUTH_BOUND], 1, exchange_field_type[Y_DIR], grid_bottom, 0,
                     &phi->val[NORTH_GUARD], 1, exchange_field_type[Y_DIR], grid_top, 0,
                     grid_comm, &status);
    }

    void simulation::exchange_charge_buffers(field *charge)
    {
        int i = 1;
        int j = 0;

        MPI_Sendrecv(&charge->val[NORTH_GUARD], 1, exchange_field_type[Y_DIR], grid_top, 0,
                     &Y_guard_data[0], N_x, MPI_DOUBLE, grid_bottom, 0, grid_comm, &status);
        if (grid_bottom != MPI_PROC_NULL)
            charge->reduceSouthBound(Y_guard_data);

        MPI_Sendrecv(&charge->val[SOUTH_GUARD], 1, exchange_field_type[Y_DIR], grid_bottom, 0,
                     &Y_guard_data[0], N_x, MPI_DOUBLE, grid_top, 0, grid_comm, &status);
        if (grid_top != MPI_PROC_NULL)
            charge->reduceNorthBound(Y_guard_data);

        MPI_Sendrecv(&charge->val[WEST_GUARD], 1, exchange_field_type[X_DIR], grid_left, 0,
                     &X_guard_data[0], N_int_y, MPI_DOUBLE, grid_right, 0, grid_comm, &status);
        if (grid_right != MPI_PROC_NULL)
            charge->reduceEastBound(X_guard_data);

        MPI_Sendrecv(&charge->val[EAST_GUARD], 1, exchange_field_type[X_DIR], grid_right, 0,
                     &X_guard_data[0], N_int_y, MPI_DOUBLE, grid_left, 0, grid_comm, &status);
        if (grid_left != MPI_PROC_NULL)
            charge->reduceWestBound(X_guard_data);
    }

    void simulation::exchange_particles_buffers(species *lepton)
    {
        // Communication to determine the size of the arrays of each buffer
        MPI_Sendrecv(&(lepton->size_send_north), 1, MPI_INT, grid_top, 0, &(lepton->size_recv_south), 1, MPI_INT, grid_bottom, 0, grid_comm, &status);
        MPI_Sendrecv(&(lepton->size_send_south), 1, MPI_INT, grid_bottom, 0, &(lepton->size_recv_north), 1, MPI_INT, grid_top, 0, grid_comm, &status);
        MPI_Sendrecv(&(lepton->size_send_west), 1, MPI_INT, grid_left, 0, &(lepton->size_recv_east), 1, MPI_INT, grid_right, 0, grid_comm, &status);
        MPI_Sendrecv(&(lepton->size_send_east), 1, MPI_INT, grid_right, 0, &(lepton->size_recv_west), 1, MPI_INT, grid_left, 0, grid_comm, &status);

        MPI_Sendrecv(&(lepton->size_send_ne), 1, MPI_INT, grid_ne, 0, &(lepton->size_recv_sw), 1, MPI_INT, grid_sw, 0, grid_comm, &status);
        MPI_Sendrecv(&(lepton->size_send_sw), 1, MPI_INT, grid_sw, 0, &(lepton->size_recv_ne), 1, MPI_INT, grid_ne, 0, grid_comm, &status);
        MPI_Sendrecv(&(lepton->size_send_nw), 1, MPI_INT, grid_nw, 0, &(lepton->size_recv_se), 1, MPI_INT, grid_se, 0, grid_comm, &status);
        MPI_Sendrecv(&(lepton->size_send_se), 1, MPI_INT, grid_se, 0, &(lepton->size_recv_nw), 1, MPI_INT, grid_nw, 0, grid_comm, &status);
        //
        // allocate memory for the vectors that are going to receive the MPI particles
        part recv_dummy;
        recv_dummy.ix = -1; // set to -1 as a way to check later if there was "actual" communication
        recv_dummy.iy = -1;
        lepton->recv_buffer_east.assign(lepton->size_recv_east, recv_dummy);
        lepton->recv_buffer_west.assign(lepton->size_recv_west, recv_dummy);
        lepton->recv_buffer_north.assign(lepton->size_recv_north, recv_dummy);
        lepton->recv_buffer_south.assign(lepton->size_recv_south, recv_dummy);

        lepton->recv_buffer_ne.assign(lepton->size_recv_ne, recv_dummy);
        lepton->recv_buffer_nw.assign(lepton->size_recv_nw, recv_dummy);
        lepton->recv_buffer_se.assign(lepton->size_recv_se, recv_dummy);
        lepton->recv_buffer_sw.assign(lepton->size_recv_sw, recv_dummy);

        //! Buffers Communication
        // All traffic in direction "top"
        MPI_Sendrecv(&(lepton->send_buffer_north[0]), lepton->send_buffer_north.size(), exchange_part_type, grid_top, 0, &(lepton->recv_buffer_south[0]), lepton->size_recv_south, exchange_part_type, grid_bottom, 0, grid_comm, &status);
        // All traf\fic in direction "bottom"
        MPI_Sendrecv(&(lepton->send_buffer_south[0]), lepton->send_buffer_south.size(), exchange_part_type, grid_bottom, 0, &(lepton->recv_buffer_north[0]), lepton->size_recv_north, exchange_part_type, grid_top, 0, grid_comm, &status);
        // All traf\fic in direction "bottom"
        MPI_Sendrecv(&(lepton->send_buffer_west[0]), lepton->send_buffer_west.size(), exchange_part_type, grid_left, 0, &(lepton->recv_buffer_east[0]), lepton->size_recv_east, exchange_part_type, grid_right, 0, grid_comm, &status);
        // All traffic in direction "right"
        MPI_Sendrecv(&(lepton->send_buffer_east[0]), lepton->send_buffer_east.size(), exchange_part_type, grid_right, 0, &(lepton->recv_buffer_west[0]), lepton->size_recv_west, exchange_part_type, grid_left, 0, grid_comm, &status);

        // All traffic in direction "ne-sw"
        MPI_Sendrecv(&(lepton->send_buffer_ne[0]), lepton->send_buffer_ne.size(), exchange_part_type, grid_ne, 0, &(lepton->recv_buffer_sw[0]), lepton->size_recv_sw, exchange_part_type, grid_sw, 0, grid_comm, &status);
        // All traf\fic in direction "sw-ne"
        MPI_Sendrecv(&(lepton->send_buffer_sw[0]), lepton->send_buffer_sw.size(), exchange_part_type, grid_sw, 0, &(lepton->recv_buffer_ne[0]), lepton->size_recv_ne, exchange_part_type, grid_ne, 0, grid_comm, &status);
        // All traf\fic in direction "se-nw"
        MPI_Sendrecv(&(lepton->send_buffer_se[0]), lepton->send_buffer_se.size(), exchange_part_type, grid_se, 0, &(lepton->recv_buffer_nw[0]), lepton->size_recv_nw, exchange_part_type, grid_nw, 0, grid_comm, &status);
        // All traffic in direction "nw-se"
        MPI_Sendrecv(&(lepton->send_buffer_nw[0]), lepton->send_buffer_nw.size(), exchange_part_type, grid_nw, 0, &(lepton->recv_buffer_se[0]), lepton->size_recv_se, exchange_part_type, grid_se, 0, grid_comm, &status);
    }

    // Jacobi solver
    void simulation::jacobi(field *phi, field *charge)
    {
        double res, e;
        double global_res = 1.0;
        double tol = 1e-7;

        long int loop = 0;

        phi->setValue(0.0);
        // Defining a new temporary field (temp is not part of the domain)
        field temp(N_x, N_y);

        if(grid_rank==4){
            phi->print_field(std::cout);
            charge->print_field(std::cout);
        }

        // Starting the iteration loop
        while (global_res > tol)
        {
            // making res 0 so that any error greater than 0 can be equated to this
            res = 0.0;

            // Making the temp field zero after every iteration
            temp.setValue(0.0);

            // exchanges buffer cells
            exchange_phi_buffers(phi);

            for (int i = 1; i <= N_int_y; i++)
                for (int j = 1; j <= N_int_x; j++)
                {
                    temp.val[POSITION] = .25 * (phi->val[NORTH] + phi->val[SOUTH] +
                                                phi->val[EAST] + phi->val[WEST] -
                                                charge->val[POSITION]);

                    e = fabs(temp.val[POSITION] - phi->val[POSITION]);
                    if (e > res) // norm infty: supremo
                        res = e;
                }

            // Transferring values from temp to u
            for (int i = 1; i <= N_int_y; i++)
                for (int j = 1; j <= N_int_x; j++)
                    phi->val[POSITION] = temp.val[POSITION];

            if (loop % 10 == 0) // balance to be found...
                MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MAX, grid_comm);

            loop++;
        }

        exchange_phi_buffers(phi);

         std::cout << "Maximum residual: " << res << "  | Number of iterations: " << loop << " | rank: " << grid_rank << std::endl;
    }

    void simulation::set_E_value(field *phi, field *Ex_field, field *Ey_field)
    {
        for (int i = 1; i <= N_int_y; i++)
            for (int j = 1; j <= N_int_x; j++)
            {
                Ex_field->val[POSITION] = (phi->val[WEST] - phi->val[EAST]) / (2.f * dx);
                Ey_field->val[POSITION] = (phi->val[SOUTH] - phi->val[NORTH]) / (2.f * dy);
            }

        if (grid_left == MPI_PROC_NULL)
        {
            if (bc[X_DIR] == CONDUCTIVE)
            {
                for (int i = 0; i < N_y; i++)
                {
                    Ey_field->val[WEST_GUARD] = 0;
                    Ex_field->val[WEST_GUARD] = -phi->val[WEST_BOUND] / dx;
                }
            }
        }

        if (grid_right == MPI_PROC_NULL)
        {
            if (bc[X_DIR] == CONDUCTIVE)
            {
                for (int i = 0; i < N_y; i++)
                {
                    Ey_field->val[EAST_GUARD] = 0;
                    Ex_field->val[EAST_GUARD] = phi->val[EAST_BOUND] / dx;
                }
            }
        }

        if (grid_top == MPI_PROC_NULL)
        {
            if (bc[Y_DIR] == CONDUCTIVE)
            {
                for (int j = 0; j < N_x; j++)
                {
                    Ex_field->val[NORTH_GUARD] = 0;
                    Ey_field->val[NORTH_GUARD] = phi->val[NORTH_BOUND] / dy;
                }
            }
        }

        if (grid_bottom == MPI_PROC_NULL)
        {
            if (bc[Y_DIR] == CONDUCTIVE)
            {
                for (int j = 0; j < N_x; j++)
                {
                    Ex_field->val[SOUTH_GUARD] = 0;
                    Ey_field->val[SOUTH_GUARD] = -phi->val[SOUTH_BOUND] / dy;
                }
            }
        }

        exchange_phi_buffers(Ex_field);
        exchange_phi_buffers(Ey_field);
    }

    /*

    // bit dumb: use bool to choose u or charge
    void simulation::write_output_u(domain &subdomain, int rank, int time)
    {
        int l, m;
        int N_local = subdomain.u->N;
        int N_local_x = subdomain.u->Nx;
        int N_local_y = subdomain.u->Ny;
        std::fstream Output_file;

        std::string filename;

        // offset of cell numbering for the subdomain
        offset[X_DIR] = grid_coord[X_DIR] * (N_local_x - 2);
        offset[Y_DIR] = grid_coord[Y_DIR] * (N_local_y - 2);

        filename = "../results/subdomain_" + std::to_string(rank) + "__t_" + std::to_string(time) + ".txt";
        std::string space = "        ";
        Output_file.open(filename, std::ios::out);
        Output_file << "x position" << space << "y position" << space << "field value" << std::endl;

        for (int i = 0; i < N_local; i++)
        {
            if (subdomain.u->bc[i] == NONE)
            {
                l = i % N_local_x;
                m = (int)i / N_local_x;

                double value_x = (l + offset[X_DIR]) * (subdomain.charge->h);
                double value_y = (m + offset[Y_DIR]) * (subdomain.charge->h);

                Output_file << "    " << value_x << space << value_y << space << subdomain.u->val[i] << std::endl;
            }
        }

        Output_file.close();
    }
    void simulation::write_output_charge(domain &subdomain, int rank, int time)
    {
        std::fstream Output_file;
        std::string filename;
        int l, m;
        int N_local = subdomain.u->N;
        int N_local_x = subdomain.u->Nx;
        int N_local_y = subdomain.u->Ny;

        filename = "../results/charge_" + std::to_string(rank) + "__t_" + std::to_string(time) + ".txt";
        std::string space = "   ";

        Output_file.open(filename, std::ios::out);
        int precision = 4;

        for (int i = 0; i < N_local; i++)
        {
            if (i % (N_local_x) == 0)
                Output_file << std::endl;

            Output_file << std::setw(precision) << subdomain.charge->val[i] << space;
            // if (i % 4 == 0)
        }

        Output_file.close();
    }
    */
}