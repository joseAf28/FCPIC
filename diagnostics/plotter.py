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


####!Variables to change
results_path = "/home/jose/Desktop/FCPIC/results/"

nb_ranks = 4
<<<<<<< HEAD
counter = 700
=======
counter = 50
>>>>>>> a400d65807475516a0dc3c4337f9839133295f59


##number and orientation of each process

######################


##Read Data from results folder

#!Charge_field case
field_type = "charge_field"
charge_field = []

for i in range(0, nb_ranks):
    field_dummy = []
    readField_rank(i, field_dummy, field_type)
    charge_field.append(field_dummy)

#!Ex_field case
field_type = "charge_field"
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


# print(Ex_field[0])
# print(len(charge_field[0]))

###
####test dimensions
# print("xlen: " )
# print(len(Ex_field[0][1]))
# print("ylen:")
# print(len(Ex_field[0][1][0]))


##animation
snapshots_charge = []
snapshots_Ex = []
snapshots_Ey = []
<<<<<<< HEAD

image_counter = 0
=======
>>>>>>> a400d65807475516a0dc3c4337f9839133295f59
for count_plot in range(0, counter):
    big_charge_dummy = []
    Ex_dummy = []
    Ey_dummy = []
    
    for i in range(0, len(charge_field[2][count_plot])):
        charge_aux = charge_field[2][count_plot][i] + charge_field[3][count_plot][i]
        big_charge_dummy.append(charge_aux)

        Ex_aux = Ex_field[2][count_plot][i] + Ex_field[3][count_plot][i]
        Ex_dummy.append(Ex_aux)

        Ey_aux = Ey_field[2][count_plot][i] + Ey_field[3][count_plot][i]
        Ey_dummy.append(Ey_aux)


<<<<<<< HEAD
    # for i in range(len(charge_field[0][count_plot])-1, 0, -1):
    for i in range(0, len(charge_field[0][count_plot])):
=======
    for i in range(len(charge_field[0][count_plot])-1, 0, -1):
>>>>>>> a400d65807475516a0dc3c4337f9839133295f59
        charge_aux = charge_field[0][count_plot][i] + charge_field[1][count_plot][i]
        big_charge_dummy.append(charge_aux)

        Ex_aux = Ex_field[0][count_plot][i] + Ex_field[1][count_plot][i]
        Ex_dummy.append(Ex_aux)

        Ey_aux = Ey_field[0][count_plot][i] + Ey_field[1][count_plot][i]
        Ey_dummy.append(Ey_aux)


    # print("*********")
    # print(big_charge[0])
    # print("*********")
    # print(big_charge[-1]) 

    # print(big_charge[-1])
    snapshots_charge.append(big_charge_dummy)
    snapshots_Ex.append(Ex_dummy)
    snapshots_Ey.append(Ey_dummy)

<<<<<<< HEAD

    # plt.figure(image_counter)
    # plt.imshow(big_charge_dummy, interpolation ='nearest')
    # plt.xlabel(r"$x$")
    # plt.ylabel(r"$y$")
    # plt.title(r"$\rho$")
    # plt.colorbar()
    # plt.savefig(results_path + "plots/charge_field_2species_"+ str(count_plot) +"_2.png")

    # image_counter = image_counter + 1

    # plt.figure(image_counter)
    # plt.imshow(Ex_dummy, interpolation ='nearest')
    # plt.xlabel(r"$x$")
    # plt.ylabel(r"$y$")
    # plt.title(r"$E_x$")
    # plt.colorbar()
    # plt.savefig(results_path + "plots/Ex_field_2species_"+ str(count_plot) +"_2.png")

    # image_counter = image_counter + 1

    # plt.figure(image_counter)
    # plt.imshow(Ey_dummy, interpolation ='nearest')
    # plt.xlabel(r"$x$")
    # plt.ylabel(r"$y$")
    # plt.title(r"$E_y$")
    # plt.colorbar()
    # plt.savefig(results_path + "plots/Ey_field_2species_"+ str(count_plot) +"_2.png")

    # image_counter = image_counter + 1

##Video Creation
fps = 30
nSeconds = math.floor(counter/fps)
print(nSeconds)
=======
    # plt.imshow(big_charge,  cmap="Greens", interpolation ='nearest')
    # plt.xlabel(r"$x$")
    # plt.ylabel(r"$y$")
    # plt.title(r"$\rho$")
    
    # plt.savefig(results_path + "/plots/test_"+ str(count_plot) +"_2.png")

fps = 10
nSeconds = math.floor(counter/fps)
print(nSeconds)
# snapshots = [ np.random.rand(5,5) for _ in range( nSeconds * fps ) ]

>>>>>>> a400d65807475516a0dc3c4337f9839133295f59
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure( figsize=(8,8) )

a = snapshots_charge[0]
<<<<<<< HEAD
im = plt.imshow(a, interpolation='spline16')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$\rho$")
plt.colorbar()
=======
im = plt.imshow(a, interpolation='nearest', cmap="Greens")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$\rho$")
>>>>>>> a400d65807475516a0dc3c4337f9839133295f59

def animate_func(i):
    if i % fps == 0:
        print( '.', end ='' )

    im.set_array(snapshots_charge[i])
    return [im]

anim = animation.FuncAnimation(fig, animate_func, 
                               frames = nSeconds * fps,
                               interval = 1000 / fps, # in ms
                               )

<<<<<<< HEAD
anim.save('charge_2species_1ppc_700_mod_vert_square_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

print('Charge Anim Done!')
=======
anim.save('test_charge_anim.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

print('Charge Anim Done!')









##Good example I found online

# import numpy as np
# def func(x, y):
#     return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

# rng = np.random.default_rng()
# points = rng.random((1000, 2))
# values = func(points[:,0], points[:,1])

# grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

# from scipy.interpolate import griddata
# grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
# grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
# grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')

# import matplotlib.pyplot as plt
# plt.subplot(221)
# plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
# plt.plot(points[:,0], points[:,1], 'k.', ms=1)
# plt.title('Original')
# plt.subplot(222)
# plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
# plt.title('Nearest')
# plt.subplot(223)
# plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
# plt.title('Linear')
# plt.subplot(224)
# plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
# plt.title('Cubic')
# plt.gcf().set_size_inches(6, 6)
# plt.show()
>>>>>>> a400d65807475516a0dc3c4337f9839133295f59
