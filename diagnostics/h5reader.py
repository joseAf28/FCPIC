import h5py

def readH5(filename, Ex_field, Ey_field, charge_field, vec_particles):
    f = h5py.File(filename, "r")

    first_key = list(f.keys())
    Ex_key = list(f[first_key[0]].keys())
    Ey_key = list(f[first_key[1]].keys())
    charge_key = list(f[first_key[2]].keys())
    particles_key = list(f[first_key[3]].keys())

    Ex_field_aux = []
    Ey_field_aux = []
    charge_field_aux = []
    vec_particles_aux = []

    for i in range(0, len(Ex_key)):
        Ex_field_aux.append(f[first_key[0]][Ex_key[i]][()])

    for i in range(0, len(Ey_key)):
        Ey_field_aux.append(f[first_key[1]][Ey_key[i]][()])

    for i in range(0, len(Ex_key)):
        charge_field_aux.append(f[first_key[2]][charge_key[i]][()])

    for i in range(0, len(Ex_key)):
        vec_particles_aux.append(f[first_key[3]][particles_key[i]][()])

    Ex_field.append(Ex_field_aux)
    Ey_field.append(Ey_field_aux)
    charge_field.append(charge_field_aux)
    vec_particles.append(vec_particles_aux)


##! varibles to change 
results_path = "/home/jose/Desktop/FCPIC/results/"
number_ranks = 6
##! 

filename_vec = []

Ex_field = []
Ey_field = []
charge_field = []
vec_particles = []

for i in range(0, number_ranks):
    filename = results_path + "data_rank_" + str(i) + ".h5"
    filename_vec.append(filename)

print(filename_vec)

for i in range(0, number_ranks):
    readH5(filename_vec[i], Ex_field, Ey_field, charge_field, vec_particles)



##debugging  hdf5 results
print(len(Ex_field[0]))
print("**************")
print(len(vec_particles))
print("**************")
print(len(vec_particles[3][1]))
print("**************")
# print(Ex_field[-1][0])
print("**************")
print(Ex_field[0][1])