import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 

def readH5(filename, Ex_field, Ey_field, charge_field, rank_id, vec_part):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    Ex_key = list(f[first_key[0]].keys())
    Ey_key = list(f[first_key[1]].keys())
    charge_key = list(f[first_key[2]].keys())
    particles_key = []
    for i in range(3, len(first_key)-1):
        particles_key_aux = list(f[first_key[i]].keys())
        particles_key.append(particles_key_aux)

    rank_id_key = list(f[first_key[-1]].keys())
    Ex_field_aux = []
    Ey_field_aux = []
    charge_field_aux = []
    vec_part_aux = []

    for i in range(0, len(Ex_key)):
        Ex_field_aux.append(f[first_key[0]][Ex_key[i]][()])

    for i in range(0, len(Ey_key)):
        Ey_field_aux.append(f[first_key[1]][Ey_key[i]][()])

    for i in range(0, len(charge_key)):
        charge_field_aux.append(f[first_key[2]][charge_key[i]][()])

    for i in range(0, len(particles_key)):
        veci_particles_aux = []
        for j in range(0, len(particles_key[i])):
                veci_particles_aux.append(f[first_key[i+3]][particles_key[i][j]][()])
        vec_part_aux.append(veci_particles_aux)


    Ex_field.append(Ex_field_aux)
    Ey_field.append(Ey_field_aux)
    charge_field.append(charge_field_aux)
    rank_id.append(f[first_key[4]][rank_id_key[0]][()])
    vec_part.append(vec_part_aux)


def H5readEx(filename, Ex_field, counter):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    Ex_key = list(f[first_key[0]].keys())
    Ex_field.append(f[first_key[0]][Ex_key[counter]][()])

def H5readEy(filename, Ey_field, counter):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    Ey_key = list(f[first_key[1]].keys())
    Ey_field.append(f[first_key[1]][Ey_key[counter]][()])

def H5readCharge(filename, Ey_field, counter):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    Ey_key = list(f[first_key[2]].keys())
    Ey_field.append(f[first_key[2]][Ey_key[counter]][()])

def H5readParticles(filename, vec_part, counter):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    particles_key = []

    for i in range(3, len(first_key)-1):
        particles_key_aux = list(f[first_key[i]].keys())
        particles_key.append(particles_key_aux)

    vec_part_aux = []
    for i in range(0, len(particles_key)):
        veci_particles_aux = (f[first_key[i+3]][particles_key[i][counter]][()])
        vec_part_aux.append(veci_particles_aux)

    vec_part.append(vec_part_aux)

def H5readRank(filename, rank_id):
    f = h5py.File(filename, "r")
    first_key = list(f.keys())
    rank_id_key = list(f[first_key[-1]].keys())
    #rank_id.append(f[first_key[4]][rank_id_key[0]][()])
    rank_id.append(f[first_key[-1]][rank_id_key[0]][()])


##########! varibles to change 
results_path = "../results/"
number_ranks = 4

grid_x_max = 4
grid_y_max = 1

counter = 60

#lx = 34./3.
#ly = 14./3.

lx = 45./3.
ly = 18./3.

dx = 1/3.
dy = 1/3.

bc = 1
#name_output = "electron_anim_"
name_output = "electron_anim_small_"

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

# print(filename_vec)

for i in range(0, len(filename_vec)):
    H5readRank(filename_vec[i], rank_id)
    

def Particle2Anim(file_vec, counter, x_data, y_data, vx_data, vy_data):
    vec_part = []
    for i in range(0, len(file_vec)):
        H5readParticles(file_vec[i], vec_part, counter)
    
    if bc == 1:
        for n_proc in range(0, len(vec_part)):
            for n_spec in range(0, len(vec_part[n_proc])):
                for i in range(0, len(vec_part[n_proc][n_spec])):
                    x_data.append(np.fmod(lx*(rank_id[n_proc][0]) + vec_part[n_proc][n_spec][i][0]*dx + vec_part[n_proc][n_spec][i][2],lx*grid_x_max))
                    y_data.append(np.fmod(ly*(rank_id[n_proc][1]) + vec_part[n_proc][n_spec][i][1]*dy + vec_part[n_proc][n_spec][i][3],ly*grid_y_max))
                    vx_data.append(vec_part[n_proc][n_spec][i][4])
                    vy_data.append(vec_part[n_proc][n_spec][i][5])
    if bc == 2:
        for n_proc in range(0, len(vec_part)):
            for n_spec in range(0, len(vec_part[n_proc])):
                for i in range(0, len(vec_part[n_proc][n_spec])):
                    x_data.append(lx*(rank_id[n_proc][0]) + vec_part[n_proc][n_spec][i][0]*dx + vec_part[n_proc][n_spec][i][2])
                    y_data.append(ly*(rank_id[n_proc][1]) + vec_part[n_proc][n_spec][i][1]*dy + vec_part[n_proc][n_spec][i][3])
                    vx_data.append(vec_part[n_proc][n_spec][i][4])
                    vy_data.append(vec_part[n_proc][n_spec][i][5])

                # print(vec_part[n_proc][n_spec][i][0])


list_indx = []
for indy in range(grid_y_max-1, -1, -1):
    list_indx_x = []
    for indx in range(0, grid_x_max):
        for n_proc in range(0, grid_x_max*grid_y_max):
            if (rank_id[n_proc][1] == indy and rank_id[n_proc][0] == indx):
                list_indx_x.append(n_proc)
    
    list_indx.append(list_indx_x)

print(rank_id)


def Ex2Anim(file_vec, counter, Ex_data):
    Ex_field = []
    for i in range(0, len(file_vec)):
        H5readEx(file_vec[i], Ex_field, counter)

    for indy in range(0, grid_y_max):
        for i in range(len(Ex_field[2])-2, 0, -1):
            Ex_aux = np.concatenate((Ex_field[list_indx[indy][0]][i][1:-2], Ex_field[list_indx[indy][1]][i][1:-2]))
            Ex_data.append(Ex_aux)


def Ey2Anim(file_vec, counter, Ey_data):
    Ey_field = []
    for i in range(0, len(file_vec)):
        H5readEy(file_vec[i], Ey_field, counter)

    for indy in range(0, grid_y_max):
        for i in range(len(Ey_field[2])-2, 0, -1):
            Ey_aux = np.concatenate((Ey_field[list_indx[indy][0]][i][1:-2], Ey_field[list_indx[indy][1]][i][1:-2]))
            Ey_data.append(Ey_aux)


def Charge2Anim(file_vec, counter, charge_data):
    charge_field = []
    for i in range(0, len(file_vec)):
        H5readCharge(file_vec[i], charge_field, counter)

    for indy in range(0, grid_y_max):
        for i in range(len(charge_field[2])-2, 0, -1):
            # charge_aux = np.concatenate((charge_field[list_indx[indy][0]][i][1:-2], charge_field[list_indx[indy][1]][i][1:-2]))

            charge_aux = (charge_field[list_indx[indy][0]][i][1:-2])
            charge_data.append(charge_aux)


def animate_particles(counter):
    x_data = []
    y_data = []
    vx_data = []
    vy_data = []
    Particle2Anim(filename_vec, counter, x_data, y_data, vx_data, vy_data)
    im.set_offsets(np.transpose([x_data, y_data]))
    return [im]


def animate_xphase(counter):
    x_data = []
    y_data = []
    vx_data = []
    vy_data = []
    Particle2Anim(filename_vec, counter, x_data, y_data, vx_data, vy_data)
    im.set_offsets(np.transpose([x_data, vx_data]))
    return [im]


def animate_yphase(counter):
    x_data = []
    y_data = []
    vx_data = []
    vy_data = []
    Particle2Anim(filename_vec, counter, x_data, y_data, vx_data, vy_data)
    im.set_offsets(np.transpose([y_data, vy_data]))
    return [im]


def animate_Exfield(counter):
    Ex_data = []
    Ex2Anim(filename_vec, counter, Ex_data)
    im.set_array(Ex_data)
    return [im]


def animate_Eyfield(counter):
    Ey_data = []
    Ey2Anim(filename_vec, counter, Ey_data)
    im.set_array(Ey_data)
    return [im]


def animate_charge(counter):
    charge_data = []
    Charge2Anim(filename_vec, counter, charge_data)
    im.set_array(charge_data)
    return [im]

snapshot_x = []
snapshot_y = []
snapshot_vx = []
snapshot_vy = []
Particle2Anim(filename_vec, 0, snapshot_x, snapshot_y, snapshot_vx, snapshot_vy)

charge_data = []
Charge2Anim(filename_vec, 1, charge_data)



# # ##! Particle animation
fps = 20
nSeconds = math.floor(counter/fps)
print(nSeconds)

fig = plt.figure()
im = plt.scatter(snapshot_x, snapshot_y, marker=".")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$Particles$")
anim = animation.FuncAnimation(fig, animate_particles, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )

anim.save(results_path+ "videos/" + name_output + "part_clean.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])

print('Particles Anim Done!')

# #X PHASE SPACE

fig = plt.figure( figsize=(8,8) )
snapshot_x = []
snapshot_y = []
snapshot_vx = []
snapshot_vy = []
Particle2Anim(filename_vec, 0, snapshot_x, snapshot_y, snapshot_vx, snapshot_vy)
im = plt.scatter(snapshot_x, snapshot_vx, marker=".")
plt.xlabel(r"$x$")
plt.ylabel(r"$vx$")
plt.title(r"$X Phase Space$")
anim = animation.FuncAnimation(fig, animate_xphase, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )
anim.save(results_path+ "videos/" + name_output + "xphase_clean.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])

print('X Phase Space Anim Done!')

#Y PHASE SPACE
fig = plt.figure( figsize=(8,8) )
im = plt.scatter(snapshot_y, snapshot_vy, marker=".")
plt.xlabel(r"$y$")
plt.ylabel(r"$vy$")
plt.title(r"$Y Phase Space$")
anim = animation.FuncAnimation(fig, animate_yphase, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )
anim.save(results_path+ "videos/" + name_output + "yphase_clean.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])

print('Y Phase Space Anim Done!')


# # # ##! Plots 
#     # plt.figure(count_plot)
#     # plt.imshow(big_charge_dummy, interpolation ='nearest')
#     # plt.xlabel(r"$x$")
#     # plt.ylabel(r"$y$")
#     # plt.title(r"$\rho$")
#     # plt.colorbar()
#     # plt.savefig(results_path + "plots/charge_field_2species_"+ str(count_plot) +"_2.png")

# # #     image_counter = image_counter + 1

# # #     plt.figure(image_counter)
# # #     plt.imshow(Ex_dummy, interpolation ='nearest')
# # #     plt.xlabel(r"$x$")
# # #     plt.ylabel(r"$y$")
# # #     plt.title(r"$E_x$")
# # #     plt.colorbar()
# # #     plt.savefig(results_path + "plots/Ex_field_2species_"+ str(count_plot) +"_2.png")

# # #     image_counter = image_counter + 1

# # #     plt.figure(image_counter)
# # #     plt.imshow(Ey_dummy, interpolation ='nearest')
# # #     plt.xlabel(r"$x$")
# # #     plt.ylabel(r"$y$")
# # #     plt.title(r"$E_y$")
# # #     plt.colorbar()
# # #     plt.savefig(results_path + "plots/Ey_field_2species_"+ str(count_plot) +"_2.png")

# # # # ##!Charge Animation
fps = 10
nSeconds = math.floor(counter/fps)
print(nSeconds)

fig = plt.figure( figsize=(6,8) )

charge_data = []
Charge2Anim(filename_vec, 0, charge_data)
a = charge_data
im = plt.imshow(a, interpolation='spline16')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$\rho$")
plt.colorbar()
anim = animation.FuncAnimation(fig, animate_charge, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )
anim.save(results_path+"videos/" + name_output + "charge_clean.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])
print('Charge Anim Done!')

# ##!Ex_field Animation
# fps = 10
# nSeconds = math.floor(counter/fps)
# print(nSeconds)

# fig = plt.figure( figsize=(8,8) )
# a = []
# Ex2Anim(filename_vec, 0, a)
# im = plt.imshow(a, interpolation='spline16')
# plt.xlabel(r"$x$")
# plt.ylabel(r"$y$")
# plt.title(r"$E_x$")

# anim = animation.FuncAnimation(fig, animate_Exfield, 
#                                frames = nSeconds * fps,
#                                interval = 1000 / fps, # in ms
#                                )
# anim.save(results_path+"videos/" + name_output + "Ex_field_clean.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])
# print('Ex field Anim Done!')

# # ##!Ey_field Animation
# fps = 10
# nSeconds = math.floor(counter/fps)
# print(nSeconds)

# a = []
# Ey2Anim(filename_vec, 0, a)
# im = plt.imshow(a, interpolation='spline16')
# plt.xlabel(r"$x$")
# plt.ylabel(r"$y$")
# plt.title(r"$E_y$")
# plt.colorbar()
# anim = animation.FuncAnimation(fig, animate_Eyfield, 
#                                frames = nSeconds * fps,
#                                interval = 1000 / fps, # in ms
#                                )
# anim.save(results_path+"videos/" + name_output + "Ey_field_clean.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])
# print('Ey field Anim Done!')


