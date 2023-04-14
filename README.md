# FCPIC

The **FCPIC** is a *2D Electrostatic Particle In Cell code with a 2D MPI parallelization*. It is written in **C++**.

To run the code, it is required to have the *hdf5* package,  the *make*, and the *g++* compiler installed.
 To run the results's analysis scripts, it is required to have *python3* installed

## Folder Structure

The code has the following folder structure:
folder src where all *.cpp* and *.h* files with the code are defined;
folder main where the objects are instantiated
- folder **lib** to where the makefile sends the formed library
- folder **bin** to where the makefile sends the executables created;
- folder **diagnostics**, where  there are *python* Scripts to analyze the results from the simulation ;
- folder **results** to where the .h5 of the results are send;

## How to compile and run the executables

To run the code, it is necessary to run the following commands in the terminal the FCPIC folder (same as the makefile):

``make filename.exe``

where filename is the name of .cpp file that instantiates the classes. This file should be in the main folder
run the following command in the bin folder

``mpirun ./filename.exe`` 

There is already a cpp file, called main.cpp in the mai folder with everything needed to compile and run the program.


## Initial conditions set up

To change the  initial conditions for the simulation, there is an input file called test.txt where there are examples of already given initial conditions:

i

1. -npart=nn1,nn2,...    Particle number for all species
2. -rand=rr1,rr2,...     Type of distribution (0->uniform grid, 1->uniform random)
3. -charge=qq1,qq2,...   Charge of all species (in q_proton/m_electron)
4. -mass=mm1,mm2,...     Mass of all species (in m_electron)
5. -temp=tt1,tt2,...     The temperature of all species (in eV)
6. -vxfluid=vx1,vx2,...  X component of fluid velocity (in v_thermal)
7. -vyfluid=vy1,vy2,...  Y component of fluid velocity (in v_thermal)
8. -nxlen=lll            Horizontal length of the simulation box (in m)
9. -nxproc=ppp           Number of MPI processes horizontally in the grid
10. -aspect=aaa           Box aspect ratio
11. -simtime=ttt          Simulation time (in secs)
12. -boundcond=bbb        Boundary condition (1->periodic, 2->conductive)

To access the help menu of the code to consult all the possible setup configurations, we need to run the executable command in the bin folder with the flag ``-help`` at the end of the command