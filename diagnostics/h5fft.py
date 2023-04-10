import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def readH5(filename, Ex_field, Ey_field, charge_field, vec_part1, vec_part2):
    f = h5py.File(filename, "r")

    first_key = list(f.keys())
    Ex_key = list(f[first_key[0]].keys())
    Ey_key = list(f[first_key[1]].keys())
    charge_key = list(f[first_key[2]].keys())

    particles_1_key = list(f[first_key[3]].keys())
    # particles_2_key = list(f[first_key[4]].keys())

    Ex_field_aux = []
    Ey_field_aux = []
    charge_field_aux = []
    vec_particles_1_aux = []
    vec_particles_2_aux = []

    for i in range(0, len(Ex_key)):
        Ex_field_aux.append(f[first_key[0]][Ex_key[i]][()])

    for i in range(0, len(Ey_key)):
        Ey_field_aux.append(f[first_key[1]][Ey_key[i]][()])

    for i in range(0, len(Ex_key)):
        charge_field_aux.append(f[first_key[2]][charge_key[i]][()])

    for i in range(0, len(Ex_key)):
        vec_particles_1_aux.append(f[first_key[3]][particles_1_key[i]][()])
    
    # for i in range(0, len(Ex_key)):
    #     vec_particles_2_aux.append(f[first_key[4]][particles_2_key[i]][()])

    Ex_field.append(Ex_field_aux)
    Ey_field.append(Ey_field_aux)
    charge_field.append(charge_field_aux)
    vec_part1.append(vec_particles_1_aux)
    # vec_part2.append(vec_particles_2_aux)


##########! varibles to change 
results_path = "/home/jose/Desktop/FCPIC/results/"
number_ranks = 4
counter = 100

counter_space = 700
lx = 16./3.
ly = 16./3.

dx = 1/3
dy = 1/3

##########! 

filename_vec = []

Ex_field = []
Ey_field = []
charge_field = []
vec_part1 = []
vec_part2 = []

for i in range(0, number_ranks):
    filename = results_path + "newdata_rank_" + str(i) + ".h5"
    filename_vec.append(filename)

# print(filename_vec)

for i in range(0, number_ranks):
    readH5(filename_vec[i], Ex_field, Ey_field, charge_field, vec_part1, vec_part2)

snapshots_charge = []
snapshots_Ex = []
snapshots_Ey = []

### rank - counter - line_field - cells that count
# image_counter = 0
for count_plot in range(0, len(Ex_field[0])):
    big_charge_dummy = []
    Ex_dummy = []
    Ey_dummy = []
    phase_dummy = []

    for i in range(len(charge_field[2][count_plot])-1, 0, -1):
    # for i in range(0, len(charge_field[2][count_plot])):
        charge_aux = np.concatenate((charge_field[2][count_plot][i][0:-2], charge_field[3][count_plot][i][0:-2]))
        big_charge_dummy.append(charge_aux)

        Ex_aux = np.concatenate((Ex_field[2][count_plot][i][0:-2], Ex_field[3][count_plot][i][0:-2]))
        Ex_dummy.append(Ex_aux)

        Ey_aux = np.concatenate((Ey_field[2][count_plot][i][0:-2], Ey_field[3][count_plot][i][0:-2]))
        Ey_dummy.append(Ey_aux)

    for i in range(len(charge_field[0][count_plot])-1, 0, -1):
    # for i in range(0, len(charge_field[0][count_plot])):
        charge_aux = np.concatenate((charge_field[0][count_plot][i][0:-2], charge_field[1][count_plot][i][0:-2]))
        big_charge_dummy.append(charge_aux)

        Ex_aux = np.concatenate((Ex_field[0][count_plot][i][0:-2], Ex_field[1][count_plot][i][0:-2]))
        Ex_dummy.append(Ex_aux)

        Ey_aux = np.concatenate((Ey_field[0][count_plot][i][0:-2], Ey_field[1][count_plot][i][0:-2]))
        Ey_dummy.append(Ey_aux)

    snapshots_charge.append(big_charge_dummy)
    snapshots_Ex.append(Ex_dummy)
    snapshots_Ey.append(Ey_dummy)

# [i][j][k] ; i - time, j - y cell, k - x cell

print(len(snapshots_Ex))
print("***************************")
print(len(snapshots_Ex[1]))
print("***************************")
print(len(snapshots_Ex[0][1]))
print("***************************")
print(len(snapshots_Ex[0][10]))
print("***************************")
print(len(snapshots_Ex))


xfield_Ex_wind = []

for i in range(0, len(snapshots_Ex)):
    xfield_Ex_wind.append(snapshots_Ex[i][10])

# print(xfield_Ex_wind[0])
# print(xfield_Ex_wind[10])

plt.plot(xfield_Ex_wind[0])
plt.plot(xfield_Ex_wind[10])
plt.plot(xfield_Ex_wind[20])
plt.plot(xfield_Ex_wind[30])
plt.savefig("ex_test.png")

fftSol = np.abs(np.fft.fft2(xfield_Ex_wind))

N = 5e-6
T =.0000000001

dspatial = N/len(snapshots_Ex[10][1])
dtime =  T/len(snapshots_Ex)
spatialfrequencies = np.fft.fftfreq(len(snapshots_Ex[10][1]), d=dspatial)
timefrequencies = np.fft.fftfreq(len(snapshots_Ex), d=dtime)


kmax = spatialfrequencies.max()
wmax = timefrequencies.max()

print("kmax")
print(kmax)

print("wmax")
print(wmax*1e-9)



# print(dspatial)
# print(1./dspatial)

# print("************")

# print(dtime)
# print(1./dtime)

# print("omega")
# print(timefrequencies)

index_vec = round(len(fftSol[1, :])/2)
c = plt.imshow(fftSol[::-1, :-index_vec], origin="lower", extent=[0, 10, 0, 10], vmin = 0, vmax = 1000)
plt.colorbar(c)
plt.xlabel(r'$k $')
plt.ylabel(r'$\omega$')
plt.savefig("structure2_S.png")
