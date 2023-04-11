import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 

def readH5(filename, Ex_field, Ey_field, charge_field, rank_id, vec_part):
    f = h5py.File(filename, "r")

    first_key = list(f.keys())
    print(first_key)

    Ex_key = list(f[first_key[0]].keys())
    Ey_key = list(f[first_key[1]].keys())
    charge_key = list(f[first_key[2]].keys())

    particles_key = []

    for i in range(3, len(first_key)-1):
        particles_key_aux = list(f[first_key[i]].keys())
        particles_key.append(particles_key_aux)


    # print(particles_key[0])
    # print(len(particles_key))

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
    # vec_part2.append(vec_particles_2_aux)


##########! varibles to change 
results_path = "../results/"
number_ranks = 4
counter = 500
counter_space = 100
lx = 16./3.
ly = 16./3.

dx = 1/3
dy = 1/3

bc = 1
name_output = "electron_anim_"

grid_x_max = 2
grid_y_max = 2
##########! 

filename_vec = []
Ex_field = []
Ey_field = []
charge_field = []

rank_id = []
vec_part = []

for i in range(0, number_ranks):
    filename = results_path + "final_sim_rank_" + str(i) + ".h5"
    filename_vec.append(filename)

# print(filename_vec)

for i in range(0, number_ranks):
    readH5(filename_vec[i], Ex_field, Ey_field, charge_field, rank_id, vec_part)


print(rank_id)
snapshot_x = []
snapshot_y = []
snapshot_vx = []
snapshot_vy = []


## vec_part index
### [i][j][k][l]] --  i: process; j: species; k: time step; l: index_part; m: particles' feature 

## rank_id index
### [i][j] -- i : process; j ::  0 - y direction;  1 - x direction 


for count_plot in range(0, len(vec_part[0][0])):
    x_data = []
    y_data = []
    vx_data = []
    vy_data = []
    
    for n_proc in range(0, len(vec_part)):
        for n_spec in range(0, len(vec_part[n_proc])):
            for i in range(0, len(vec_part[n_proc][n_spec][count_plot])):
                x_data.append(lx*(rank_id[n_proc][1]) + vec_part[n_proc][n_spec][count_plot][i][0]*dx + vec_part[n_proc][n_spec][count_plot][i][2])
                y_data.append(ly*(rank_id[n_proc][0]) + vec_part[n_proc][n_spec][count_plot][i][1]*dy + vec_part[n_proc][n_spec][count_plot][i][3])
                vx_data.append(vec_part[n_proc][n_spec][count_plot][i][4])
                vy_data.append(vec_part[n_proc][n_spec][count_plot][i][5])


    snapshot_x.append(x_data)
    snapshot_y.append(y_data)
    snapshot_vx.append(vx_data)
    snapshot_vy.append(vy_data)

# # ##! Particle animation
fps = 20
nSeconds = math.floor(len(vec_part[0][0])/fps)
print(nSeconds)
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure( figsize=(8,8) )

im = plt.scatter(snapshot_x[0], snapshot_y[0], marker=".")


plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$Particles$")

def animate_func(i):
    if i % fps == 0:
        print( '.', end ='' )

    im.set_offsets(np.transpose([snapshot_x[i], snapshot_y[i]]))
    return [im]

anim = animation.FuncAnimation(fig, animate_func, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )

anim.save(results_path+ "videos/" + name_output + "part.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])

print('Particles Anim Done!')

#X PHASE SPACE

fig = plt.figure( figsize=(8,8) )

im = plt.scatter(snapshot_x[0], snapshot_vx[0], marker=".")

plt.xlabel(r"$x$")
plt.ylabel(r"$vx$")
plt.title(r"$X Phase Space$")

def animate_func(i):
    if i % fps == 0:
        print( '.', end ='' )

    im.set_offsets(np.transpose([snapshot_x[i], snapshot_vx[i]]))
    return [im]

anim = animation.FuncAnimation(fig, animate_func, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )

anim.save(results_path+ "videos/" + name_output + "xphase.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])

print('X Phase Space Anim Done!')

#Y PHASE SPACE

fig = plt.figure( figsize=(8,8) )

im = plt.scatter(snapshot_y[0], snapshot_vy[0], marker=".")

plt.xlabel(r"$y$")
plt.ylabel(r"$vy$")
plt.title(r"$Y Phase Space$")

def animate_func(i):
    if i % fps == 0:
        print( '.', end ='' )

    im.set_offsets(np.transpose([snapshot_y[i], snapshot_vy[i]]))
    return [im]

anim = animation.FuncAnimation(fig, animate_func, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )

anim.save(results_path+ "videos/" + name_output + "yphase.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])

print('Y Phase Space Anim Done!')



# ##!! FIELDS CASE
snapshots_charge = []
snapshots_Ex = []
snapshots_Ey = []


##order index

list_indx = []
for indy in range(grid_y_max-1, -1, -1):
    list_indx_x = []
    for indx in range(0, grid_x_max):
        for n_proc in range(0, 4):
            if (rank_id[n_proc][0] == indy and rank_id[n_proc][1] == indx):
                list_indx_x.append(n_proc)
    
    list_indx.append(list_indx_x)

print(list_indx)


### rank - counter - line_field - cells that count
# image_counter = 0
for count_plot in range(0, len(charge_field[2])):
    big_charge_dummy = []
    Ex_dummy = []
    Ey_dummy = []
    phase_dummy = []

    ##!!!!manual concatenatenation in the right grid structure in x direction
    for indy in range(0, grid_y_max):
        for i in range(len(charge_field[2][count_plot])-2, 0, -1):
            charge_aux = np.concatenate((charge_field[list_indx[indy][0]][count_plot][i][1:-2], charge_field[list_indx[indy][1]][count_plot][i][1:-2]))
            big_charge_dummy.append(charge_aux)

            Ex_aux = np.concatenate((Ex_field[list_indx[indy][0]][count_plot][i][1:-2], Ex_field[list_indx[indy][1]][count_plot][i][1:-2]))
            Ex_dummy.append(Ex_aux)

            Ey_aux = np.concatenate((Ey_field[list_indx[indy][0]][count_plot][i][1:-2], Ey_field[list_indx[indy][1]][count_plot][i][1:-2]))
            Ey_dummy.append(Ey_aux)


    snapshots_charge.append(big_charge_dummy)
    snapshots_Ex.append(Ex_dummy)
    snapshots_Ey.append(Ey_dummy)
# # ##! Plots 
    # plt.figure(count_plot)
    # plt.imshow(big_charge_dummy, interpolation ='nearest')
    # plt.xlabel(r"$x$")
    # plt.ylabel(r"$y$")
    # plt.title(r"$\rho$")
    # plt.colorbar()
    # plt.savefig(results_path + "plots/charge_field_2species_"+ str(count_plot) +"_2.png")

# #     image_counter = image_counter + 1

# #     plt.figure(image_counter)
# #     plt.imshow(Ex_dummy, interpolation ='nearest')
# #     plt.xlabel(r"$x$")
# #     plt.ylabel(r"$y$")
# #     plt.title(r"$E_x$")
# #     plt.colorbar()
# #     plt.savefig(results_path + "plots/Ex_field_2species_"+ str(count_plot) +"_2.png")

# #     image_counter = image_counter + 1

# #     plt.figure(image_counter)
# #     plt.imshow(Ey_dummy, interpolation ='nearest')
# #     plt.xlabel(r"$x$")
# #     plt.ylabel(r"$y$")
# #     plt.title(r"$E_y$")
# #     plt.colorbar()
# #     plt.savefig(results_path + "plots/Ey_field_2species_"+ str(count_plot) +"_2.png")

# #     image_counter = image_counter + 1

# # ##!Charge Animation
fps = 10
nSeconds = math.floor(len(charge_field[2])/fps)
print(nSeconds)
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure( figsize=(6,8) )

a = snapshots_charge[0]
im = plt.imshow(a, interpolation='spline16')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$\rho$")
plt.colorbar()
def animate_func(i):
    if i % fps == 0:
        print( '.', end ='' )
    im.set_array(snapshots_charge[i])
    return [im]
anim = animation.FuncAnimation(fig, animate_func, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )
anim.save(results_path+"videos/" + name_output + "charge.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])
print('Charge Anim Done!')

##!Ex_field Animation
fps = 10
nSeconds = math.floor(len(charge_field[2])/fps)
print(nSeconds)
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure( figsize=(8,8) )
a = snapshots_Ex[0]
im = plt.imshow(a, interpolation='spline16')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$E_x$")
def animate_func(i):
    if i % fps == 0:
        print( '.', end ='' )
    im.set_array(snapshots_Ex[i])
    return [im]
anim = animation.FuncAnimation(fig, animate_func, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )
anim.save(results_path+"videos/" + name_output + "Ex_field.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])
print('Ex field Anim Done!')

##!Ey_field Animation
fps = 10
nSeconds = math.floor(len(charge_field[2])/fps)
print(nSeconds)
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure( figsize=(8,8) )
a = snapshots_Ey[0]
im = plt.imshow(a, interpolation='spline16')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$E_y$")
plt.colorbar()
def animate_func(i):
    if i % fps == 0:
        print( '.', end ='' )
    im.set_array(snapshots_Ey[i])
    return [im]
anim = animation.FuncAnimation(fig, animate_func, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )
anim.save(results_path+"videos/" + name_output + "Ey_field.mp4", fps=fps, extra_args=['-vcodec', 'libx264'])
print('Ey field Anim Done!')


