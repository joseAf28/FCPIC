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
counter = 40
counter_space = 400
lx = 21
ly = 11

dx = 1.;
dy = 1.;
##########! 

filename_vec = []

Ex_field = []
Ey_field = []
charge_field = []
vec_part1 = []
vec_part2 = []

for i in range(0, number_ranks):
    filename = results_path + "data_rank_" + str(i) + ".h5"
    filename_vec.append(filename)

# print(filename_vec)

for i in range(0, number_ranks):
    readH5(filename_vec[i], Ex_field, Ey_field, charge_field, vec_part1, vec_part2)

# snapshot_x = []
# snapshot_y = []

# for count_plot in range(0, counter_space):
#     x_data = []
#     y_data = []

#     for i in range(0, len(vec_part1[0][count_plot])):
#         x_data.append(vec_part1[0][count_plot][i][0]*dx + vec_part1[0][count_plot][i][2])
#         y_data.append(vec_part1[0][count_plot][i][1]*dy + vec_part1[0][count_plot][i][3])

#     for i in range(0, len(vec_part1[1][count_plot])):
#         x_data.append(lx - 1 + vec_part1[1][count_plot][i][0]*dx + vec_part1[1][count_plot][i][2])
#         y_data.append(vec_part1[1][count_plot][i][1]*dy + vec_part1[1][count_plot][i][3])

#     for i in range(0, len(vec_part1[2][count_plot])):
#         x_data.append(vec_part1[2][count_plot][i][0]*dx + vec_part1[2][count_plot][i][2])
#         y_data.append(ly - 1 + vec_part1[2][count_plot][i][1]*dy + vec_part1[2][count_plot][i][3])

#     for i in range(0, len(vec_part1[3][count_plot])):
#         x_data.append(lx -1 + vec_part1[3][count_plot][i][0]*dx + vec_part1[3][count_plot][i][2])
#         y_data.append(ly -1 + vec_part1[3][count_plot][i][1]*dy + vec_part1[3][count_plot][i][3])


#     for i in range(0, len(vec_part2[0][count_plot])):
#         x_data.append(vec_part2[0][count_plot][i][0]*dx + vec_part2[0][count_plot][i][2])
#         y_data.append(vec_part2[0][count_plot][i][1]*dy + vec_part2[0][count_plot][i][3])

#     for i in range(0, len(vec_part2[1][count_plot])):
#         x_data.append(lx + vec_part2[1][count_plot][i][0]*dx + vec_part2[1][count_plot][i][2])
#         y_data.append(vec_part2[1][count_plot][i][1]*dy + vec_part2[1][count_plot][i][3])

#     for i in range(0, len(vec_part2[2][count_plot])):
#         x_data.append(vec_part2[2][count_plot][i][0]*dx + vec_part2[2][count_plot][i][2])
#         y_data.append(ly + vec_part2[2][count_plot][i][1]*dy + vec_part2[2][count_plot][i][3])

#     for i in range(0, len(vec_part2[3][count_plot])):
#         x_data.append(lx + vec_part2[3][count_plot][i][0]*dx + vec_part2[3][count_plot][i][2])
#         y_data.append(ly + vec_part2[3][count_plot][i][1]*dy + vec_part2[3][count_plot][i][3])


#     snapshot_x.append(x_data)
#     snapshot_y.append(y_data)


##!! FIELDS CASE
# ##animation
snapshots_charge = []
snapshots_Ex = []
snapshots_Ey = []

### rank - counter - line_field - cells that count
# image_counter = 0
for count_plot in range(0, counter):
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

# print(len(snapshots_charge[0]))
# print("***************************")
# print(len(snapshots_charge[1]))
# print("***************************")
# print(len(snapshots_charge[0][1]))
# print("***************************")
# print(len(snapshots_charge[0][10]))
# print("***************************")
# print(len(snapshots_charge))

# plt.plot(snapshots_charge[0][1])
# plt.plot(snapshots_charge[200][1])
# plt.savefig("test.png")

xfield_charge_wind = []

len_wind = len(snapshots_charge[0][10]) 

for i in range(0, len(snapshots_charge)):

    xfield_charge_wind.append(snapshots_charge[i][10]*np.hanning(len_wind))


N = 40
T = 40

fftSol = np.abs(np.fft.fft2(xfield_charge_wind))
tpCountX = len(xfield_charge_wind[1])
tpCount = T
dspatial = len(xfield_charge_wind[1])/N
dtime = len(xfield_charge_wind[:][1])/T

spatialfrequencies = np.fft.fftfreq(tpCountX, d=dspatial)
timefrequencies = np.fft.fftfreq(tpCount, d=dtime)

kmax = spatialfrequencies.max()
wmax = timefrequencies.max()


# print(len(fftSol))
# print("***********************")
# print(len(fftSol[0]))
# print("***********************")
# print(len(fftSol[0][0]))
# print("***********************")
# print(fftSol[0][0][0])


size_array = len(fftSol[1, :])
index_vec = round(size_array/2)

# c = plt.imshow(fftSol, extent =[0, kmax, 0, wmax],origin ='lower', vmin=0, interpolation="spline16")


c = plt.imshow(fftSol[::-1, :-index_vec], extent =[0, kmax, 0, wmax],origin ='lower', vmin=0, vmax=1000)

k_vec = np.linspace(0, wmax, 30)
w_vec = np.linspace(0, wmax, 30)

plt.colorbar(c)
plt.xlabel(r'$k $')
plt.ylabel(r'$\omega$')
plt.savefig("structure2_S.png")

print(fftSol[1,:-index_vec])
print("**********************")
print(fftSol[1])

print("**********************")
print(fftSol[:,1])

