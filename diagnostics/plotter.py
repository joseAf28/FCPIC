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
counter = 100
counter_space = 700
lx = 21
ly = 11

##number and orientation of each process

######################

###!! SPACE CASE

#!Phase Space case
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


print(len(phase_0_field[0][0][:]))


print(len(phase_1_field[0][:]))

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



## Particle animation
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

anim.save(results_path+'videos/particles_anim_4_off.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

print('Particles Anim Done!')



###!! FIELDS CASE

# #!Charge_field case
# field_type = "charge_field"
# charge_field = []

# for i in range(0, nb_ranks):
#     field_dummy = []
#     readField_rank(i, field_dummy, field_type)
#     charge_field.append(field_dummy)

# #!Ex_field case
# field_type = "Ex_field"
# Ex_field = []

# for i in range(0, nb_ranks):
#     field_dummy = []
#     readField_rank(i, field_dummy, field_type)
#     Ex_field.append(field_dummy)

# #!Ey_field case
# field_type = "Ey_field"
# Ey_field = []

# for i in range(0, nb_ranks):
#     field_dummy = []
#     readField_rank(i, field_dummy, field_type)
#     Ey_field.append(field_dummy)



# print(snapshots_phase[1])
# print(snapshots_phase[-1])

# print(Ex_field[0])
# print(len(charge_field[0]))

# plt.figure()
# plt.scatter(snapshot_x[0], snapshot_y[0], marker=".")
# plt.savefig("test_part_0.png")

# plt.figure()
# plt.scatter(snapshot_x[-1], snapshot_y[-1], marker=".")
# plt.savefig("test_part_-1.png")


###
####test dimensions
# print("xlen: " )
# print(len(Ex_field[0][1]))
# print("ylen:")
# print(len(Ex_field[0][1][0]))


# # ##animation
# snapshots_charge = []
# snapshots_Ex = []
# snapshots_Ey = []
# snapshots_Phase = []

# # image_counter = 0
# for count_plot in range(0, counter):
#     big_charge_dummy = []
#     Ex_dummy = []
#     Ey_dummy = []
#     phase_dummy = []
    
#     for i in range(0, len(charge_field[2][count_plot])):
#         charge_aux = charge_field[2][count_plot][i] + charge_field[3][count_plot][i]
#         big_charge_dummy.append(charge_aux)

#         Ex_aux = Ex_field[2][count_plot][i] + Ex_field[3][count_plot][i]
#         Ex_dummy.append(Ex_aux)

#         Ey_aux = Ey_field[2][count_plot][i] + Ey_field[3][count_plot][i]
#         Ey_dummy.append(Ey_aux)


#     # for i in range(len(charge_field[0][count_plot])-1, 0, -1):
#     for i in range(0, len(charge_field[0][count_plot])):
#         charge_aux = charge_field[0][count_plot][i] + charge_field[1][count_plot][i]
#         big_charge_dummy.append(charge_aux)

#         Ex_aux = Ex_field[0][count_plot][i] + Ex_field[1][count_plot][i]
#         Ex_dummy.append(Ex_aux)

#         Ey_aux = Ey_field[0][count_plot][i] + Ey_field[1][count_plot][i]
#         Ey_dummy.append(Ey_aux)

#     for i in range(0, len(charge_field[2][count_plot])):

    


#     # print("*********")
#     # print(big_charge[0])
#     # print("*********")
#     # print(big_charge[-1]) 

#     # print(big_charge[-1])
#     snapshots_charge.append(big_charge_dummy)
#     snapshots_Ex.append(Ex_dummy)
#     snapshots_Ey.append(Ey_dummy)


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

# #Video Creation
# fps = 30
# nSeconds = math.floor(counter/fps)
# print(nSeconds)
# # First set up the figure, the axis, and the plot element we want to animate
# fig = plt.figure( figsize=(8,8) )

# a = snapshots_charge[0]
# im = plt.imshow(a, interpolation='spline16')
# plt.xlabel(r"$x$")
# plt.ylabel(r"$y$")
# plt.title(r"$\rho$")
# plt.colorbar()

# def animate_func(i):
#     if i % fps == 0:
#         print( '.', end ='' )

#     im.set_array(snapshots_charge[i])
#     return [im]

# anim = animation.FuncAnimation(fig, animate_func, 
#                                frames = nSeconds * fps,
#                                interval = 1000 / fps, # in ms
#                                )

# anim.save(results_path+'videos/test_charge_2species_1ppc_700_mod_vert_square_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

# print('Charge Anim Done!')
