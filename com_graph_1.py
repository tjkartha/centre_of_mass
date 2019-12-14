# ----------------------------- com_graph_1.py -------------------------------
# This file accepts a .gro file (or any single frame data in .gro format) as
# input and outputs a .csv file, which is converted to a numpy array to
# calculate centre of masses. This program assumes that you have two different
# residues A and B. The user gets to decide what A and B should be. The output
# is a print of the centre of mass arrays of A and B, along with an 
# interactive 3D plot of the same. Also, the centres of masses of each of the 
# residues will be saved as separate .txt files.

# Importing required packages
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# This function formats the .gro file into a .csv file. 
def formatter():
	my_file = raw_input("Enter the input file: ")
	os.system("cp " + my_file + " test_1.gro")
	output_file = raw_input("Enter the ouput file name with a .csv extension (other than \"test_1.csv\"): ")
	print "Chopping heads >>>>>"
	os.system("sed -i 1,2d test_1.gro")
	print "Chopping tails >>>>>"
	os.system("sed -i '$d' test_1.gro")
	print "Column separations >>>>>"
	os.system("sed -i -e 's/.\{44\}/& /'  test_1.gro")
	os.system("sed -i -e 's/.\{36\}/& /'  test_1.gro")
	os.system("sed -i -e 's/.\{28\}/& /'  test_1.gro")
	os.system("sed -i -e 's/.\{20\}/& /'  test_1.gro")
	os.system("sed -i -e 's/.\{15\}/& /'  test_1.gro")
	os.system("sed -i -e 's/.\{5\}/& /'  test_1.gro") 
	os.system("tr -s ' ' < test_1.gro | tr ' ' ',' > test_1.csv")
	os.system("mv test_1.csv " + output_file)
	os.system("rm test_1.gro")


# Mass inputs
def mass_handler(mol_frame):
	total_size = int(len(mol_frame))
	mol_num = int(mol_frame[-1, 0]) - int(mol_frame[0, 0]) + 1
	mol_size = total_size / mol_num
	atom_types = mol_frame[:mol_size, 2]
	atomic_mass_list = []
	print "Enter the atomic masses of atoms in residue ", str(mol_frame[0, 1]), " :"
	print "---------------------------------------------------------"
	for x in atom_types:
		y = float(raw_input("Mass of atom " + x + " : "))
		atomic_mass_list.append(y)
	print "---------------------------------------------------------"
	atomic_mass_array = np.array(atomic_mass_list)
	total_mass = np.sum(atomic_mass_array)
	atomic_mass_array_full = np.tile(atomic_mass_array, mol_num)
	atomic_mass_array_full = atomic_mass_array_full.reshape((len(atomic_mass_array_full), 1))
	mol_with_mass = np.hstack((mol_frame, atomic_mass_array_full))
	return mol_with_mass, mol_size, total_mass, mol_num


# Calculating the centre of masses
def mass_xyz(mol_with_mass, total_mass, num, mol_size):
	mass_x = []
	mass_y = []
	mass_z = []
	for x in mol_with_mass:
		mass_x.append(float(x[4])*float(x[7]))
		mass_y.append(float(x[5])*float(x[7]))	
		mass_z.append(float(x[6])*float(x[7]))
	mass_x = np.array(mass_x)
	#print mass_x
	mass_y = np.array(mass_y)
	mass_z = np.array(mass_z)
	mol_mass_xyz = np.hstack((mol_with_mass, mass_x.reshape((len(mass_x), 1)), mass_y.reshape((len(mass_y), 1)), mass_z.reshape((len(mass_z), 1))))
	mass_only = mol_mass_xyz[:, 7:].astype(float)
	com_list_x = []
	com_list_y = []
	com_list_z = []
	counter = 1
	while counter <= num:
		a = np.sum(mass_only[((counter-1)*mol_size):(counter*mol_size), 1]) / total_mass
		com_list_x.append(a)
		b = np.sum(mass_only[((counter-1)*mol_size):(counter*mol_size), 2]) / total_mass
		com_list_y.append(b)
		c = np.sum(mass_only[((counter-1)*mol_size):(counter*mol_size), 3]) / total_mass
		com_list_z.append(c)
		counter += 1	
	com_list_x = np.array(com_list_x)
	com_list_y = np.array(com_list_y)
	com_list_z = np.array(com_list_z)
	com_list = np.hstack((com_list_x.reshape(len(com_list_x), 1), com_list_y.reshape(len(com_list_y), 1), com_list_z.reshape(len(com_list_z), 1)))
	mol_com = np.hstack((np.sort(np.unique(mol_with_mass[:,:1]).astype(int)).reshape(len(com_list), 1), mol_with_mass[:len(com_list), 1].reshape(len(com_list), 1), com_list))
	return mol_com


formatter()

frame = np.genfromtxt("1.csv", delimiter=',', dtype=str)
frame = frame[:, 1:8]

A = raw_input("Enter the label of residue A: ")
B = raw_input("Enter the label of residue B: ")

A_frame = frame[frame[:, 1] == str(A)]
B_frame = frame[frame[:, 1] == str(B)]

A_num = 0
B_num = 0
A_with_mass, A_size, A_total_mass, A_num = mass_handler(A_frame)
B_with_mass, B_size, B_total_mass, B_num = mass_handler(B_frame)

A_com = mass_xyz(A_with_mass, A_total_mass, A_num, A_size)
B_com = mass_xyz(B_with_mass, B_total_mass, B_num, B_size)

print "\nCentres of mass of ", str(A), " :\n"
print A_com
print "\n------------------------\n"
print "Centres of mass of ", str(B), " :\n"
print B_com
print "\n------------------------\n"

np.savetxt(str(A)+".txt", A_com[:, 2:].astype(float), fmt='%.8f   %.8f   %.8f')
np.savetxt(str(B)+".txt", A_com[:, 2:].astype(float), fmt='%.8f   %.8f   %.8f')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(A_com[:, 2].astype(float), A_com[:, 3].astype(float), A_com[:, 4].astype(float), c='r', marker='o')
ax.scatter(B_com[:, 2].astype(float), B_com[:, 3].astype(float), B_com[:, 4].astype(float), c='b', marker='o')
plt.show()

print "\n-------------- Exited -----------------\n"


