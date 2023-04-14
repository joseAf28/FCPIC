//FCPIC - 2D Particle-in-Cell code using MPI
//Guilherme Crispim, João Palma, José Afonso
//Advanced Topics in Computational Physics, 2023, IST

//File simulation_io.cc:
//Implementation of all input handling, parameter setting, 
//printing and output related functions in the class Simulation

#include "simulation.hh"
#include <fstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <unistd.h>

namespace FCPIC
{
    /* 
        Function print_SI
        Inputs: numeric value, precision
        Output: string with value written using SI prefixes (e.g 10^3 -> k)
        and provided number of decimal places.
        + Used for printing in more legible manner all numeric info in the
        title screen
    */
    std::string simulation::print_SI(float x, int precision)
    {
        //Condition for negative numbers 
        //(otherwise they mess with the log10 function)
        int sign = 0;
        if (x < 0) 
        {
            sign = 1;
            x = -x;
        }

        //Separating mantissa and order of magnitude
        int exponent = x == 0 ? 0 : floor(log10(x));
        exponent -= (exponent % 3 + 3) % 3;
        float mantissa = x / pow(10., exponent);

        if (sign == 1)
            mantissa = -mantissa;

        std::string output = std::to_string(mantissa);
        output = output.substr(0, output.find(".") + precision + 1); //Precision setting

        if (exponent == -15)
            output.append(" f");
        if (exponent == -12)
            output.append(" p");
        if (exponent == -9)
            output.append(" n");
        if (exponent == -6)
            output.append(" μ");
        if (exponent == -3)
            output.append(" m");
        if (exponent == 0)
            output.append(" ");
        if (exponent == 3)
            output.append(" k");
        if (exponent == 6)
            output.append(" M");
        if (exponent == 9)
            output.append(" G");
        if (exponent == 12)
            output.append(" T");
        if (exponent == 15)
            output.append(" P");

        return output;
    }

    /* 
        Function printTitle
        + Prints a summary of all important simulation parameters
    */
    void simulation::printTitle()
    {
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
        std::cout << "Simulation size: " << print_SI(xlen * Lref, 3) << "m x " << print_SI(xlen * aspect * Lref, 3) << "m\n";
        std::cout << "Spatial discretization: " << print_SI(dx * Lref, 3) << "m x "
                  << print_SI(dy * Lref, 3) << "m"
                  << " (" << N_total_x << "x" << N_total_y << " cells)\n";
        std::cout << "Simulation time: " << print_SI(simtime * Tref, 3) << "s\n";
        std::cout << "Time discretization: " << print_SI(dt * Tref, 3) << "s ("
                  << (int)(simtime / dt) << " time steps)\n";
        std::cout << "MPI process grid: " << n_Procs << " processes ("
                  << grid[X_DIR] << "x" << grid[Y_DIR] << ")\n";
        std::cout << "Boundary condition: " << ((bc == 1) ? "periodic\n" : "conductive\n");
        std::cout << "Electron Debye length: " << print_SI(Lref, 3) << "m\n";
        std::cout << "Plasma frequency: " << print_SI(1 / Tref, 3) << "Hz\n";
        std::cout << "\n";
        for (int i = 0; i < Nspecies; i++)
        {
            std::cout << "Species " << i << ": " << Npart[i] << " particles, ";
            std::cout << ((rand_true[i] == 1) ? "random uniform dist.\n  " : "evenly distributed\n  ");
            std::cout << "q = " << charge[i] << " e, m = " << mass[i] << " me, T = " << print_SI(temp[i] * Tempref, 3) << "eV";
            std::cout << " , fluid velocity = (" << print_SI(vxfluid[i] * Vref, 3) << "m/s, " << print_SI(vyfluid[i] * Vref, 3) << "m/s)\n";
        }
        std::cout << "\n";
    }

    /* 
        Function printProgress
        Inputs: fraction of the progress of the simulation (between 0 and 1)
        + Prints a summary of all important simulation parameters
    */
    void simulation::printProgress(float prog)
    {
        int perc = round(prog * 100);
        static bool finished = false;
        if (!finished) //Test for not rewriting the progress once the sim is finished
        {
            std::cout << "\r\033[?25l";
            std::cout << "RUNNING: │";

            //Cycle writing bar length corrrespondent to progress fraction
            for (int k = 0; k < 100; k += 4)
            {
                if (perc - k >= 4)
                    std::cout << "█";
                else if (perc - k == 2 || perc - k == 3)
                    std::cout << "▌";
                else
                    std::cout << " ";
            }
            std::cout << "│ " << perc << "%"<< std::flush;
            if (perc >= 100)
            {
                std::cout << "\n\n";
                std::cout << "Simulation finished successfully";
                std::cout << "\n\n";
                std::cout << "\033[?25h";
                finished = true;
            }
        }
    }

    /* 
        Function printTime
        Inputs: file name for writing all time data
        + Prints a table summarizing the times each part of the simulation took
        (calculating fields, advancing particles, etc.)
        + Writes all this data in a separate txt file
    */
    void simulation::printTime(std::string filename)
    {
        if (grid_rank == 0)
        {
            std::cout << "  Procs  │  Sim. Setup  │ Particle Push │  Field Solve  │  HDF5 Write  ║   TOTAL";
            std::cout << std::endl;
            std::ofstream outfile(filename + ".txt");
            outfile << "# Procs: Sim. Setup (s), Particle Push (s), Field Solve (s), HDF5 Write (s), TOTAL (s)\n\n";
            outfile.close();
        }
        //Each process writes (by rank order) its section of the table and of the file
        for (int i = 0; i < n_Procs; i++)
        {
            MPI_Barrier(grid_comm);
            if (grid_rank == i)
            {
                std::cout << "─────────┼──────────────┼───────────────┼───────────────┼──────────────╫────────────\n";
                std::cout << "   " << std::right << std::setw(3) << i << "   │";
                std::cout << " T:" << std::right << std::setw(8) << setup_time * 1000 << "ms │";
                std::cout << " T:" << std::right << std::setw(10) << particle_time << "s │";
                std::cout << " T:" << std::right << std::setw(10) << field_time << "s │";
                std::cout << " T:" << std::right << std::setw(9) << hdf5_time << "s ║";
                std::cout << " " << std::right << std::setw(9) << total_time << "s";
                std::cout << std::endl;
                std::cout << std::setprecision(3);
                std::cout << "         │";
                std::cout << " %: " << std::right << std::setw(8) << setup_time * 100. / total_time << "% │";
                std::cout << " %: " << std::right << std::setw(9) << particle_time * 100. / total_time << "% │";
                std::cout << " %: " << std::right << std::setw(9) << field_time * 100. / total_time << "% │";
                std::cout << " %: " << std::right << std::setw(8) << hdf5_time * 100. / total_time << "% ║";
                std::cout << std::endl;

                std::ofstream outfile(filename + ".txt", std::ios_base::app);
                outfile << i << ": ";
                outfile << setup_time << " , ";
                outfile << particle_time << " , ";
                outfile << field_time << " , ";
                outfile << hdf5_time << " , ";
                outfile << total_time << "\n";
                outfile.close();
            }
            MPI_Barrier(grid_comm);
        }

        if(grid_rank == 0){
            std::cout << "\n";
            std::cout << "Process time report written in " << filename << ".txt" << std::endl;
        }
    }

    /* 
        Function printHelp
        + Prints a help menu explaining how to run the program and how to set all
        relevant inputs
    */
    void simulation::printHelp()
    {
        std::cout << "\n";
        std::cout << ">> mpiexec ./FCPIC.exec -infile=infile.txt\n";
        std::cout << ">> mpiexec ./FCPIC.exec -npart=1000,2000 -charge=1,-1  ...\n";
        std::cout << "\n";
        std::cout << "-infile=fff           Read parameters from input file\n";
        std::cout << "\n";
        std::cout << "-npart=nn1,nn2,...    Particle number for all species\n";
        std::cout << "                            Default: 1000\n";
        std::cout << "-rand=rr1,rr2,...     Type of distribution (0->uniform grid, 1->uniform random)\n";
        std::cout << "                            Default: 1\n";
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
        std::cout << "-aspect=aaa           Box aspect ratio (ylen = aspect * xlen)\n";
        std::cout << "                            Default: 1\n";
        std::cout << "-simtime=ttt          Simulation time (in secs)\n";
        std::cout << "                            Default: 1\n";
        std::cout << "-boundcond=bbb        Boundary condition (1->periodic, 2->conductive)\n";
        std::cout << "                            Default: 1\n";
        std::cout << "\n";
    }

    /* 
        Function readArgs
        Inputs: int main argument number and list
        + Reads the main arguments and sets them internally in the simulation
    */
    void simulation::readArgs(int argc, char **argv)
    {
        std::vector<std::string> allArgs(argv, argv + argc); //Vector with all args
        std::vector<std::string> numbers;
        std::string line, header, number;       //Aux variables for cycles
        std::vector<bool> def_values(12, true); //Flags for setting default values
        int k;

        //Cycle over each arg: numeric chars go to "number" string, letters to header.
        //If numberic chars are separated by commas, they are separated as different entries
        //in a vector
        for (auto &arg : allArgs)
        {
            k = 0;
            for (char &c : arg)
            {
                if ((c > 47 && c < 58) || c == '.' || (c == '-' && k != 0))
                    number.push_back(c);
                if (c == ',')
                {
                    numbers.push_back(number);
                    number.clear();
                }
                if ((c > 96 && c < 123) || (c == '-' && k == 0))
                    header.push_back(c);
                if (c > 64 && c < 91)
                    header.push_back(c + 32);
                
                //Only case in which no numeric values are expected next
                if (header.compare("-infile") == 0)
                {
                    arg.erase(0, 8);
                    break;
                }
                k++;
            }

            if (header.compare("-help") == 0 && rank == 0) //Flag for help menu
            {
                printHelp();
            }

            if (header.compare("-infile") == 0) //Flag for file read
            {
                def_values.assign(12, true);         //All previous values are neglected
                getParamsfromFile(arg, &def_values); //Reading file
                break; 
            }

            if (header.compare("-npart") == 0) //Flag for particle number
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    Npart.push_back(std::stoi(num));
                def_values[0] = false;
            }
            if (header.compare("-charge") == 0) //Flag for species charge
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    charge.push_back(std::stod(num));
                def_values[1] = false;
            }
            if (header.compare("-mass") == 0) //Flag for species mass
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    mass.push_back(std::stod(num));
                def_values[2] = false;
            }
            if (header.compare("-temp") == 0) //Flag for species temperature
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    temp.push_back(std::stod(num));
                def_values[3] = false;
            }
            if (header.compare("-vxfluid") == 0) //Flag for species fluid velocity (X)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    vxfluid.push_back(std::stod(num));
                def_values[4] = false;
            }
            if (header.compare("-vyfluid") == 0) //Flag for species fluid velocity (Y)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    vyfluid.push_back(std::stod(num));
                def_values[5] = false;
            }
            if (header.compare("-xlen") == 0) //Flag for simulation length in X direction
            {
                xlen = stod(number);
                def_values[6] = false;
            }
            if (header.compare("-nxproc") == 0) //Flag for number of processes along X direction
            {
                grid[X_DIR] = stoi(number);
                def_values[7] = false;
            }
            if (header.compare("-aspect") == 0) //Flag for box aspect ratio
            {
                aspect = stod(number);
                def_values[8] = false;
            }
            if (header.compare("-simtime") == 0) //Flag for simulation time
            {
                simtime = stod(number);
                def_values[9] = false;
            }
            if (header.compare("-boundcond") == 0) //Flag for boundary condition
            {
                bc = stoi(number);
                def_values[10] = false;
            }
            if (header.compare("-rand") == 0) //Flag for species distribution type
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    rand_true.push_back(std::stod(num));
                def_values[11] = false;
            }

            numbers.clear();
            number.clear();
            header.clear();
        }

        //Set default values whenever needed
        if (def_values[0])
            Npart.push_back(1000);
        if (def_values[1])
            for (auto &v : Npart)
                charge.push_back(-1.);
        if (def_values[2])
            for (auto &v : Npart)
                mass.push_back(1.);
        if (def_values[3])
            for (auto &v : Npart)
                temp.push_back(1.);
        if (def_values[4])
            for (auto &v : Npart)
                vxfluid.push_back(0.);
        if (def_values[5])
            for (auto &v : Npart)
                vyfluid.push_back(0.);
        if (def_values[6])
            xlen = 1;
        if (def_values[7])
            grid[X_DIR] = n_Procs / 2;
        if (def_values[8])
            aspect = 1;
        if (def_values[9])
            simtime = 1;
        if (def_values[10])
            bc = PERIODIC;
        if (def_values[11])
            for (auto &v : Npart)
                rand_true.push_back(1);
    }

    /* 
        Function getParamsfromFile
        Inputs: input file name, vector of flags for setting default values
        (as needed in function readArgs)
        + Reads the arguments from a file and sets them internally
    */
    void simulation::getParamsfromFile(std::string filename, std::vector<bool> *def_values)
    {

        std::ifstream infile(filename);

        if (!infile)
            {
                if (rank == 0)
                    std::cout << "ERROR OPENING INPUT FILE\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        
        std::vector<std::string> filelines, numbers;
        std::string line, header, number;

        while (getline(infile, line))
            filelines.push_back(line);

        infile.close();

        //Same cycle format as in readArgs, changing only header format
        for (auto &fline : filelines)
        {
            for (char &c : fline)
            {
                if ((c > 47 && c < 58) || c == '.' || c == '-')
                {
                    number.push_back(c);
                }
                if (c == ',')
                {
                    numbers.push_back(number);
                    number.clear();
                }
                if (c > 96 && c < 123)
                    header.push_back(c);
                if (c > 64 && c < 91)
                    header.push_back(c + 32);
                if (c == '#')
                    break;
            }

            if (header.compare("npart") == 0)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    Npart.push_back(std::stoi(num));
                (*def_values)[0] = false;
            }
            if (header.compare("charge") == 0)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    charge.push_back(std::stod(num));
                (*def_values)[1] = false;
            }
            if (header.compare("mass") == 0)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    mass.push_back(std::stod(num));
                (*def_values)[2] = false;
            }
            if (header.compare("temp") == 0)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    temp.push_back(std::stod(num));
                (*def_values)[3] = false;
            }
            if (header.compare("vxfluid") == 0)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    vxfluid.push_back(std::stod(num));
                (*def_values)[4] = false;
            }
            if (header.compare("vyfluid") == 0)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    vyfluid.push_back(std::stod(num));
                (*def_values)[5] = false;
            }
            if (header.compare("xlen") == 0)
            {
                xlen = stod(number);
                (*def_values)[6] = false;
            }
            if (header.compare("nxproc") == 0)
            {
                grid[X_DIR] = stoi(number);
                (*def_values)[7] = false;
            }
            if (header.compare("aspect") == 0)
            {
                aspect = stod(number);
                (*def_values)[8] = false;
            }
            if (header.compare("simtime") == 0)
            {
                simtime = stod(number);
                (*def_values)[9] = false;
            }
            if (header.compare("boundcond") == 0)
            {
                bc = stoi(number);
                (*def_values)[10] = false;
            }
            if (header.compare("rand") == 0)
            {
                numbers.push_back(number);
                for (auto &num : numbers)
                    rand_true.push_back(std::stod(num));
                (*def_values)[11] = false;
            }

            numbers.clear();
            number.clear();
            header.clear();
        }
    }

    /* 
        Function setParams
        + Processes the inputs given for a consistent simulation
        + Normalizes all quantities
    */
    void simulation::setParams()
    {
        //Sequence of tests checking the validity of the inputs

        if (Npart.size() != charge.size() || Npart.size() != mass.size() ||
            Npart.size() != temp.size() || Npart.size() != vxfluid.size() ||
            Npart.size() != vyfluid.size() || Npart.size() != rand_true.size())
        {
            if (rank == 0)
                std::cout << "ERROR: PLEASE PROVIDE ARGS FOR ALL SPECIES\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        Nspecies = Npart.size();

        for (int k = 0; k < Nspecies; k++)
        {
            if (Npart[k] < 1)
            {
                if (rank == 0)
                    std::cout << "ERROR: NON-PHYSICAL NUMBER OF PARTICLES (SPECIES NO " << k << ")\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            if (mass[k] <= 0)
            {
                if (rank == 0)
                    std::cout << "ERROR: NON-PHYSICAL MASS (SPECIES NO " << k << ")\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            if (temp[k] < 0)
            {
                if (rank == 0)
                    std::cout << "ERROR: NON-PHYSICAL TEMPERATURE (SPECIES NO " << k << ")\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            if (rand_true[k] != 0 && rand_true[k] != 1)
            {
                if (rank == 0)
                    std::cout << "ERROR: INVALID PARTICLE DISTRIBUTION TYPE (SPECIES NO " << k << ")\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        if (xlen <= 0)
        {
            if (rank == 0)
                std::cout << "ERROR: NON-PHYSICAL BOX LENGTH\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (n_Procs % grid[X_DIR] != 0)
        {
            if (rank == 0)
                std::cout << "ERROR: UNABLE TO FORM A PROCESS GRID\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (aspect <= 0)
        {
            if (rank == 0)
                std::cout << "ERROR: INVALID ASPECT RATIO\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (simtime <= 0)
        {
            if (rank == 0)
                std::cout << "ERROR: INVALID SIMULATION TIME\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (bc != 1 && bc != 2)
        {
            if (rank == 0)
                std::cout << "ERROR: INVALID BOUNDARY CONDITIONS\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //In case an evenly distributed grid of particles is asked,
        //the total number of particles is slightly changed in case
        //no integer distribution across X and Y can be found 
        for (int i = 0; i < Nspecies; i++)
        {
            if (rand_true[i] == 0)
            {
                //For a uniform grid, we assume Nypart = aspect * Nxpart
                Nxpart.push_back(round(sqrt(((float)Npart[i]) / aspect)));
                Nypart.push_back(round(sqrt(((float)Npart[i]) * aspect)));
                Npart[i] = Nxpart[i] * Nypart[i]; //(Possibly) changed total number of particles
            }
        }

        //Normalization constants
        Vref = 419382.88 * sqrt(temp[0]);                 //Normalizing speeds to Vthermal of 1st species (m/s)
        Nref = (float)Npart[0] / (aspect * xlen * xlen);  //Normalizing densities to 2D Density (m^-2)
        Lref = sqrt(55263494 * temp[0] / pow(Nref, 1.5)); //Normalizing lengths to electron Debye length (m)
        Tref = Lref / Vref;                               //Normalizing times to electron inverse wpe (s)
        Tempref = temp[0];                                //Normalizing temperature to 1st species T

        //Normalizing all temperatures and thermal speeds
        for (int k = 0; k < Nspecies; k++)
        {
            temp[k] /= Tempref;
            vxfluid[k] *= sqrt(temp[k]);
            vyfluid[k] *= sqrt(temp[k]);
        }

        //Normalizing times and lengths
        simtime /= Tref;
        xlen /= Lref;

        //Spatial discretization is set initially to 1/3 of Debye length (stability condition)
        dx = 1. / 3.;

        //First guesses for total number of cells across X and Y
        N_total_x = round((float)xlen / (float)dx);
        N_total_y = round(aspect * (float)N_total_x);

        grid[Y_DIR] = n_Procs / grid[X_DIR]; //number of processes across Y

        //The grid will have the same dimensions for any boundary conditions. However, if the 
        //system has fixed boundaries, our code assumes it has an additional row of physical 
        //cells across X and Y (all outward guard cell rows become physical, while, in periodic
        //BCs, the two guard cell rows represent the same space). Hence, this condition will 
        //remove 1 unit from the number of cells across both directions, so that dividing it by 
        //the number of processes should be an integer (the numberof cells per process per direction)
        if (bc != PERIODIC)
        {
            N_total_x -= 1;
            N_total_y -= 1;
        }

        //The number of cells in the X direction is adjusted slightly so that its division across 
        //all processes is evenly made
        int k = 0;
        while (N_total_x % grid[X_DIR] != 0)
        {
            k = -k - abs(k) / k;
            N_total_x += k;
        }

        N_int_x = N_total_x / grid[X_DIR]; //Number of physical cells in X per process
        N_x = N_int_x + 2;                 //Number of physical cells in X per process + 2 guard cells

        //Same procedure for Y
        k = 0;
        while (N_total_y % grid[Y_DIR] != 0)
        {
            k = -k - abs(k) / k;
            N_total_y += k;
        }

        N_int_y = N_total_y / grid[Y_DIR]; //Number of physical cells in Y per process
        N_y = N_int_y + 2;                 //Number of physical cells in Y per process + 2 guard cells

        //Recovering the extra physical cell for non periodic BCs
        if (bc != PERIODIC)
        {
            N_total_x += 1;
            N_total_y += 1;
        }

        N = N_x * N_y; //New cell number (of a single process)

        //Adjusted dx and dy for the calculated number of cells
        dx = xlen / (float)N_total_x;
        dy = aspect * xlen / (float)N_total_y;

        //CFL condition. Velocity used is the max between thermal and fluid velocities
        dt = 1 / (std::max(1., sqrt(vxfluid[0] * vxfluid[0] + vyfluid[0] * vyfluid[0])) * (1 / dx + 1 / dy));

        //Important MPI parameter. If periodic, the MPI grid can connect automatically
        //boundary processes. This variable carries that setting option 
        wrap_around[X_DIR] = bc == PERIODIC ? 1 : 0;
        wrap_around[Y_DIR] = bc == PERIODIC ? 1 : 0;
    }

    /* 
        Function confirmParams
        + Asks for a final confirmation of the parameters from part of the 
        user, after the title screen printing
        + Useful, since the parameters are slightly tweaked for the good
        functioning of the system, and this option avoids having to abort
        in "less elegant" ways the program
    */
    void simulation::confirmParams()
    {
        char type, buffer[128];

        //Only one process interacts with the user
        if (grid_rank == 0)
        {
            //Cycle until the user provides a valid answer
            do
            {
                std::cout << "Do you want to progress with these parameters [y/n]? ";
                std::cin >> buffer;
                type = buffer[0];
            } while (!std::cin.fail() && type != 'y' && type != 'n' && type != 'Y' && type != 'N');

            if (type == 'y' || type == 'Y')
            {
                sim_true = 1;
                std::cout << "\n";
            }

            if (type == 'n' || type == 'N')
                sim_true = 0;
            
            //User interaction has to be done by a single process. Info needs to travel
            //to all processes, so it should be sent in an MPI message
            for (int i = 1; i < n_Procs; i++)
                MPI_Send(&sim_true, 1, MPI_INT, i, 0, grid_comm);
        }
        
        //All other processes receive the user answer via MPI
        else
            MPI_Recv(&sim_true, 1, MPI_INT, 0, 0, grid_comm, MPI_STATUS_IGNORE);
    }
}