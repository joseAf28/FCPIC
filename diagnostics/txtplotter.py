import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 


#Generic data reader from Files
def readData(filename, Data):

    df = pd.read_csv(filename, encoding='utf-8', skiprows=0)
    df_size = len(df.index)
    i = 0
    while(i < df_size):
        line = df.iloc[i]
        x = line[0].split()
        vec = []
        for k in range(0, len(x)):
            try:
                vec.append(float(x[k]))
            except:
                pass
        
        if(len(vec)>0):
            Data.append(vec)
        i = i + 1

    return df_size

#reads all time steps files for each type of file selected
def readField_rank(rank, charge_field, field_type):
    charge_name_array = []
    for i in range(0, counter):
        name = results_path + field_type + "/rank:_" + str(rank) + "_counter_" + str(i) + ".txt"
        charge_name_array.append(name)

    for i in range(0, counter):
        charge = []
        readData(charge_name_array[i], charge)
        charge_field.append(charge)

def readSpace_rank(rank, charge_field, field_type, spec):
    charge_name_array = []
    for i in range(0, counter_space):
        name = results_path + field_type + "/particles_" + str(rank) + "_spec_" + str(spec) + "_counter_" + str(i) + ".txt"
        charge_name_array.append(name)

    for i in range(0, counter_space):
        charge = []
        readData(charge_name_array[i], charge)
        charge_field.append(charge)


####!Variables to change
results_path = "/home/jose/Desktop/FCPIC/results/"

nb_spec = 2
nb_ranks = 4
counter = 200
counter_space = 700
lx = 21
ly = 11

##number and orientation of each process

######################

# ###!! SPACE CASE

# #!Phase Space case
field_type = "phase_space"
phase_0_field = []

for i in range(0, nb_ranks):
    field_dummy = []
    readSpace_rank(i, field_dummy, field_type, 0)
    phase_0_field.append(field_dummy)


field_type = "phase_space"
phase_1_field = []

for i in range(0, nb_ranks):
    field_dummy = []
    readSpace_rank(i, field_dummy, field_type, 1)
    phase_1_field.append(field_dummy)


snapshot_x = []
snapshot_y = []

for count_plot in range(0, counter_space):
    x_data = []
    y_data = []

    for i in range(0, len(phase_0_field[0][count_plot][:])):
        x_data.append(phase_0_field[0][count_plot][i][0])
        y_data.append(phase_0_field[0][count_plot][i][1])

    for i in range(0, len(phase_0_field[1][count_plot][:])):
        x_data.append(lx + phase_0_field[1][count_plot][i][0])
        y_data.append(phase_0_field[1][count_plot][i][1])

    for i in range(0, len(phase_0_field[2][count_plot][:])):
        x_data.append(phase_0_field[2][count_plot][i][0])
        y_data.append(ly + phase_0_field[2][count_plot][i][1])

    for i in range(0, len(phase_0_field[3][count_plot][:])):
        x_data.append(lx + phase_0_field[3][count_plot][i][0])
        y_data.append(ly + phase_0_field[3][count_plot][i][1])

    for i in range(0, len(phase_1_field[0][count_plot][:])):
        x_data.append(phase_1_field[0][count_plot][i][0])
        y_data.append(phase_1_field[0][count_plot][i][1])

    for i in range(0, len(phase_1_field[1][count_plot][:])):
        x_data.append(lx + phase_1_field[1][count_plot][i][0])
        y_data.append(phase_1_field[1][count_plot][i][1])

    for i in range(0, len(phase_1_field[2][count_plot][:])):
        x_data.append(phase_1_field[2][count_plot][i][0])
        y_data.append(ly + phase_1_field[2][count_plot][i][1])

    for i in range(0, len(phase_1_field[3][count_plot][:])):
        x_data.append(lx + phase_1_field[3][count_plot][i][0])
        y_data.append(ly + phase_1_field[3][count_plot][i][1])

    snapshot_x.append(x_data)
    snapshot_y.append(y_data)

##! Particle animation
fps = 30
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

anim.save(results_path+'videos/particles_5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

print('Particles Anim Done!')


###!! FIELDS CASE

#!Charge_field case
field_type = "charge_field"
charge_field = []

for i in range(0, nb_ranks):
    field_dummy = []
    readField_rank(i, field_dummy, field_type)
    charge_field.append(field_dummy)

#!Ex_field case
field_type = "Ex_field"
Ex_field = []

for i in range(0, nb_ranks):
    field_dummy = []
    readField_rank(i, field_dummy, field_type)
    Ex_field.append(field_dummy)

#!Ey_field case
field_type = "Ey_field"
Ey_field = []

for i in range(0, nb_ranks):
    field_dummy = []
    readField_rank(i, field_dummy, field_type)
    Ey_field.append(field_dummy)


# ##animation
snapshots_charge = []
snapshots_Ex = []
snapshots_Ey = []
snapshots_Phase = []

# image_counter = 0
for count_plot in range(0, counter):
    big_charge_dummy = []
    Ex_dummy = []
    Ey_dummy = []
    phase_dummy = []

        # for i in range(len(charge_field[2][count_plot])-1, 0, -1):
    for i in range(0, len(charge_field[2][count_plot])):
        charge_aux = charge_field[2][count_plot][i] + charge_field[3][count_plot][i]
        big_charge_dummy.append(charge_aux)

        Ex_aux = Ex_field[2][count_plot][i] + Ex_field[3][count_plot][i]
        Ex_dummy.append(Ex_aux)

        Ey_aux = Ey_field[2][count_plot][i] + Ey_field[3][count_plot][i]
        Ey_dummy.append(Ey_aux)
    
    # for i in range(len(charge_field[0][count_plot])-1, 0, -1):
    for i in range(0, len(charge_field[0][count_plot])):
        charge_aux = charge_field[0][count_plot][i] + charge_field[1][count_plot][i]
        big_charge_dummy.append(charge_aux)

        Ex_aux = Ex_field[0][count_plot][i] + Ex_field[1][count_plot][i]
        Ex_dummy.append(Ex_aux)

        Ey_aux = Ey_field[0][count_plot][i] + Ey_field[1][count_plot][i]
        Ey_dummy.append(Ey_aux)

    snapshots_charge.append(big_charge_dummy)
    snapshots_Ex.append(Ex_dummy)
    snapshots_Ey.append(Ey_dummy)

##! Plots 
#     # plt.figure(image_counter)
#     # plt.imshow(big_charge_dummy, interpolation ='nearest')
#     # plt.xlabel(r"$x$")
#     # plt.ylabel(r"$y$")
#     # plt.title(r"$\rho$")
#     # plt.colorbar()
#     # plt.savefig(results_path + "plots/charge_field_2species_"+ str(count_plot) +"_2.png")

#     # image_counter = image_counter + 1

#     # plt.figure(image_counter)
#     # plt.imshow(Ex_dummy, interpolation ='nearest')
#     # plt.xlabel(r"$x$")
#     # plt.ylabel(r"$y$")
#     # plt.title(r"$E_x$")
#     # plt.colorbar()
#     # plt.savefig(results_path + "plots/Ex_field_2species_"+ str(count_plot) +"_2.png")

#     # image_counter = image_counter + 1

#     # plt.figure(image_counter)
#     # plt.imshow(Ey_dummy, interpolation ='nearest')
#     # plt.xlabel(r"$x$")
#     # plt.ylabel(r"$y$")
#     # plt.title(r"$E_y$")
#     # plt.colorbar()
#     # plt.savefig(results_path + "plots/Ey_field_2species_"+ str(count_plot) +"_2.png")

#     # image_counter = image_counter + 1

##!Charge Animation
fps = 10
nSeconds = math.floor(counter/fps)
print(nSeconds)
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure( figsize=(8,8) )

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
anim.save(results_path+'videos/charge_5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])
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
anim.save(results_path+'videos/Ex_field_5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])
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
anim.save(results_path+'videos/Ey_field_5_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])
print('Ey field Anim Done!')

    # for indy in range(0, max_y):
    #     charge_aux_y = []
    #     for n_proc2 in range(0, len(vec_part)):
            
    #         if(rank_id[n_proc2][0] == indy):
    #             print(n_proc2)
    #             charge_aux_line = []

    #             for indx in range(0, max_x):
    #                 charge_aux_x = []
    #                 for n_proc in range(0, len(vec_part)):

    #                     charge_aux_x_2 = []
    #                     if(rank_id[n_proc][1] == indx):

    #                          for i in range(len(charge_field[n_proc][count_plot])-1, 0, -1):
    #                             charge_aux_x_aux = []
    #                             charge_aux_x_aux = np.concatenate((charge_field[n_proc][count_plot][i][0:-2], charge_aux_x_aux))
    #                             charge_aux_x_2.append(charge_aux_x_aux)
                            
    #                     charge_aux_x.append(charge_aux_x_2)
    #                     # if(rank_id[n_proc][1] == indx):
    #                     #     charge_aux = np.concatenate((charge_field[n_proc][count_plot][0:-2], charge_aux))

    #             charge_aux_line.append(charge_aux_x)

    #         charge_aux_y.append(charge_aux_line)

    #     big_charge_dummy.append(charge_aux_line)
           