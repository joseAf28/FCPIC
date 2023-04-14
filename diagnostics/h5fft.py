import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def H5readRank(filename, rank_id):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    rank_id_key = list(f[first_key[-1]].keys())
    rank_id.append(f[first_key[-1]][rank_id_key[0]][()])


##########! varibles to change 
results_path = "../results/"
number_ranks = 4

grid_x_max = 4
grid_y_max = 1

counter = 238

N_x = 119
N_y = 49


lx = 15.
ly = 15.

# lx = 45./3.
# ly = 18./3.

dx = 1.
dy = 1.


dt = 0.1

bc = 1
#name_output = "electron_anim_"
name_output = "final_sim_small_rank"

##########! 

filename_vec = []
Ex_field = []
Ey_field = []
charge_field = []

rank_id = []
vec_part = []

for i in range(0, number_ranks):
    filename = results_path + "final_sim_small_rank_" + str(i) + ".h5"
    filename_vec.append(filename)


for i in range(0, len(filename_vec)):
    H5readRank(filename_vec[i], rank_id)


def H5readEy(filename, Ey_real_field, N_x, N_y, counter):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    Ey_key = list(f[first_key[1]].keys())
    Ey_field = []
    Ey_field.append(f[first_key[1]][Ey_key[counter]][()])

    Ey_real_field_aux = []
    for i in range(1, N_y-1):
        Ey_real_y_aux = []
        for j in range(0, N_x):
            Ey_real_y_aux.append(Ey_field[0][j + i*N_x])
        Ey_real_field_aux.append(Ey_real_y_aux)

    Ey_real_field.append(Ey_real_field_aux)


def H5readEx(filename, Ex_real_field, N_x, N_y, counter):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    Ex_key = list(f[first_key[0]].keys())
    Ex_field = []
    Ex_field.append(f[first_key[0]][Ex_key[counter]][()])

    Ex_real_field_aux = []
    for i in range(1, N_y-1):
        Ex_real_y_aux = []
        for j in range(0, N_x):
            Ex_real_y_aux.append(Ex_field[0][j + i*N_x])
        Ex_real_field_aux.append(Ex_real_y_aux)

    Ex_real_field.append(Ex_real_field_aux)


def H5readCharge(filename, charge_real_field, N_x, N_y, counter):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    charge_key = list(f[first_key[2]].keys())
    charge_field = []
    charge_field.append(f[first_key[2]][charge_key[counter]][()])

    charge_real_field_aux = []
    for i in range(1, N_y-1):
        charge_real_y_aux = []
        for j in range(0, N_x):
            charge_real_y_aux.append(charge_field[0][j + i*N_x])
        charge_real_field_aux.append(charge_real_y_aux)

    charge_real_field.append(charge_real_field_aux)


list_indx = []
for indx in range(0, grid_x_max):
    list_indx_x = []
    for indy in range(0, grid_y_max):
        for n_proc in range(0, grid_x_max*grid_y_max):
            if (rank_id[n_proc][1] == indy and rank_id[n_proc][0] == indx):
                list_indx_x.append(n_proc)
    
    list_indx.append(list_indx_x)

print(rank_id)
print(list_indx)



def Ex2Anim(file_vec, counter, Ex_data):
    Ex_field = []
    for i in range(0, len(file_vec)):
        H5readEx(file_vec[i], Ex_field, N_x, N_y, counter)

    for indy in range(0, grid_y_max):
        for i in range(len(Ex_field[2])-2, 0, -1):            
            #4x1
            Ex_aux = np.concatenate((Ex_field[list_indx[0][indy]][i][1:-1], Ex_field[list_indx[1][indy]][i][1:-1], Ex_field[list_indx[2][indy]][i][1:-1], Ex_field[list_indx[3][indy]][i][1:-1]))

            #2x2
            # Ex_aux = np.concatenate((Ex_field[list_indx[indy][0]][i][1:-2], Ex_field[list_indx[indy][1]][i][1:-2]))

            Ex_data.append(Ex_aux)

snapshots_Ex = []

for i in range(0, counter):
    Ex_dummy_fft =  []
    Ex2Anim(filename_vec, i, Ex_dummy_fft)
    snapshots_Ex.append(Ex_dummy_fft[10])


fftSol = np.abs(np.fft.fft2(snapshots_Ex))

N = 500e-9
T =.0000000000003

dspatial = N/len(snapshots_Ex[10])
dtime =  T/len(snapshots_Ex)
spatialfrequencies = np.fft.fftfreq(len(snapshots_Ex[10]), d=dspatial)
timefrequencies = np.fft.fftfreq(len(snapshots_Ex), d=dtime)


kmax = spatialfrequencies.max()
wmax = timefrequencies.max()


index_vec = round(len(fftSol[1, :])/2)
c = plt.imshow(fftSol[::-5, :-index_vec], origin="lower", extent=[0, kmax*1e-8, 0, wmax*1e-14], vmin=0, vmax=100)
plt.colorbar(c)
plt.xlabel(r'$k [m^{-1}]   10^{8}$ ')
plt.ylabel(r'$\omega [s^{-1}]   10^{14}$')
plt.savefig("structure2_S.png")
