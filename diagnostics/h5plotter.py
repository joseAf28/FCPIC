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
    particles_2_key = list(f[first_key[4]].keys())

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
    
    for i in range(0, len(Ex_key)):
        vec_particles_2_aux.append(f[first_key[4]][particles_2_key[i]][()])

    Ex_field.append(Ex_field_aux)
    Ey_field.append(Ey_field_aux)
    charge_field.append(charge_field_aux)
    vec_part1.append(vec_particles_1_aux)
    vec_part2.append(vec_particles_2_aux)


##########! varibles to change 
results_path = "../results/"
number_ranks = 4
counter = 500
counter_space = 500
lx = 16./3.
ly = 16./3.

dx = 1/3
dy = 1/3

bc = 1
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

snapshot_x = []
snapshot_y = []
snapshot_vx = []
snapshot_vy = []

for count_plot in range(0, counter_space):
    x_data = []
    y_data = []
    vx_data = []
    vy_data = []
    
    for i in range(0, len(vec_part1[0][count_plot])):
        x_data.append(vec_part1[0][count_plot][i][0]*dx + vec_part1[0][count_plot][i][2])
        y_data.append(vec_part1[0][count_plot][i][1]*dy + vec_part1[0][count_plot][i][3])
        vx_data.append(vec_part1[0][count_plot][i][4])
        vy_data.append(vec_part1[0][count_plot][i][5])

    for i in range(0, len(vec_part1[1][count_plot])):
        x_data.append(np.fmod(lx + vec_part1[1][count_plot][i][0]*dx + vec_part1[1][count_plot][i][2],2*lx))
        y_data.append(vec_part1[1][count_plot][i][1]*dy + vec_part1[1][count_plot][i][3])
        vx_data.append(vec_part1[1][count_plot][i][4])
        vy_data.append(vec_part1[1][count_plot][i][5])

    for i in range(0, len(vec_part1[2][count_plot])):
        x_data.append(vec_part1[2][count_plot][i][0]*dx + vec_part1[2][count_plot][i][2])
        y_data.append(np.fmod(ly + vec_part1[2][count_plot][i][1]*dy + vec_part1[2][count_plot][i][3],2*ly))
        vx_data.append(vec_part1[2][count_plot][i][4])
        vy_data.append(vec_part1[2][count_plot][i][5])

    for i in range(0, len(vec_part1[3][count_plot])):
        x_data.append(np.fmod(lx + vec_part1[3][count_plot][i][0]*dx + vec_part1[3][count_plot][i][2],2*lx))
        y_data.append(np.fmod(ly + vec_part1[3][count_plot][i][1]*dy + vec_part1[3][count_plot][i][3],2*ly))
        vx_data.append(vec_part1[3][count_plot][i][4])
        vy_data.append(vec_part1[3][count_plot][i][5])
    

    for i in range(0, len(vec_part2[0][count_plot])):
        x_data.append(vec_part2[0][count_plot][i][0]*dx + vec_part2[0][count_plot][i][2])
        y_data.append(vec_part2[0][count_plot][i][1]*dy + vec_part2[0][count_plot][i][3])
        vx_data.append(vec_part2[0][count_plot][i][4])
        vy_data.append(vec_part2[0][count_plot][i][5])

    for i in range(0, len(vec_part2[1][count_plot])):
        x_data.append(lx + vec_part2[1][count_plot][i][0]*dx + vec_part2[1][count_plot][i][2])
        y_data.append(vec_part2[1][count_plot][i][1]*dy + vec_part2[1][count_plot][i][3])
        vx_data.append(vec_part2[1][count_plot][i][4])
        vy_data.append(vec_part2[1][count_plot][i][5])

    for i in range(0, len(vec_part2[2][count_plot])):
        x_data.append(vec_part2[2][count_plot][i][0]*dx + vec_part2[2][count_plot][i][2])
        y_data.append(ly + vec_part2[2][count_plot][i][1]*dy + vec_part2[2][count_plot][i][3])
        vx_data.append(vec_part2[2][count_plot][i][4])
        vy_data.append(vec_part2[2][count_plot][i][5])

    for i in range(0, len(vec_part2[3][count_plot])):
        x_data.append(lx + vec_part2[3][count_plot][i][0]*dx + vec_part2[3][count_plot][i][2])
        y_data.append(ly + vec_part2[3][count_plot][i][1]*dy + vec_part2[3][count_plot][i][3])
        vx_data.append(vec_part2[3][count_plot][i][4])
        vy_data.append(vec_part2[3][count_plot][i][5])


    snapshot_x.append(x_data)
    snapshot_y.append(y_data)
    snapshot_vx.append(vx_data)
    snapshot_vy.append(vy_data)

##! Particle animation
fps = 20
nSeconds = math.floor(counter_space/fps)
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

anim.save(results_path+'videos/newtest2_part_hdf5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

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

anim.save(results_path+'videos/newtest2_xphase_hdf5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

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

anim.save(results_path+'videos/newtest2_yphase_hdf5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

print('Y Phase Space Anim Done!')

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

# ##! Plots 
#     plt.figure(image_counter)
#     plt.imshow(big_charge_dummy, interpolation ='nearest')
#     plt.xlabel(r"$x$")
#     plt.ylabel(r"$y$")
#     plt.title(r"$\rho$")
#     plt.colorbar()
#     plt.savefig(results_path + "plots/charge_field_2species_"+ str(count_plot) +"_2.png")

#     image_counter = image_counter + 1

#     plt.figure(image_counter)
#     plt.imshow(Ex_dummy, interpolation ='nearest')
#     plt.xlabel(r"$x$")
#     plt.ylabel(r"$y$")
#     plt.title(r"$E_x$")
#     plt.colorbar()
#     plt.savefig(results_path + "plots/Ex_field_2species_"+ str(count_plot) +"_2.png")

#     image_counter = image_counter + 1

#     plt.figure(image_counter)
#     plt.imshow(Ey_dummy, interpolation ='nearest')
#     plt.xlabel(r"$x$")
#     plt.ylabel(r"$y$")
#     plt.title(r"$E_y$")
#     plt.colorbar()
#     plt.savefig(results_path + "plots/Ey_field_2species_"+ str(count_plot) +"_2.png")

#     image_counter = image_counter + 1

# ##!Charge Animation
fps = 10
nSeconds = math.floor(counter/fps)
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
anim.save(results_path+'videos/newtest2_charge_hdf5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])
print('Charge Anim Done!')

##!Ex_field Animation
fps = 10
nSeconds = math.floor(counter/fps)
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
anim.save(results_path+'videos/newtest2_Ex_field_hdf5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])
print('Ex field Anim Done!')

##!Ey_field Animation
fps = 10
nSeconds = math.floor(counter/fps)
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
anim.save(results_path+'videos/newtest2_Ey_field_hdf5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])
print('Ey field Anim Done!')


