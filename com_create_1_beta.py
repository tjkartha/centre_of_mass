import numpy as np


def mass_handler(mol_frame):
	total_size = int(len(mol_frame))
	mol_num = int(mol_frame[-1, 0]) - int(mol_frame[0, 0]) + 1
	mol_size = total_size / mol_num
	atom_types = mol_frame[:mol_size, 2]
	atomic_mass_list = []
	for x in atom_types:
		y = float(raw_input("Enter the mass of atom " + x + " : "))
		atomic_mass_list.append(y)
	atomic_mass_array = np.array(atomic_mass_list)
	total_mass = np.sum(atomic_mass_array)
	atomic_mass_array_full = np.tile(atomic_mass_array, mol_num)
	atomic_mass_array_full = atomic_mass_array_full.reshape((len(atomic_mass_array_full), 1))
	print mol_frame.shape
	print atomic_mass_array_full.shape
	mol_with_mass = np.hstack((mol_frame, atomic_mass_array_full))
	return mol_with_mass, mol_size, total_mass, mol_num


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
		#print a
		b = np.sum(mass_only[((counter-1)*mol_size):(counter*mol_size), 2]) / total_mass
		com_list_y.append(b)
		#print b
		c = np.sum(mass_only[((counter-1)*mol_size):(counter*mol_size), 3]) / total_mass
		com_list_z.append(c)
		#print c
		counter += 1	
	com_list_x = np.array(com_list_x)
	com_list_y = np.array(com_list_y)
	com_list_z = np.array(com_list_z)
	com_list = np.hstack((com_list_x.reshape(len(com_list_x), 1), com_list_y.reshape(len(com_list_y), 1), com_list_z.reshape(len(com_list_z), 1)))
	mol_com = np.hstack((np.sort(np.unique(mol_with_mass[:,:1]).astype(int)).reshape(len(com_list), 1), mol_with_mass[:len(com_list), 1].reshape(len(com_list), 1), com_list))
	return mol_com


frame = np.genfromtxt("1.csv", delimiter=',', dtype=str)

frame = frame[:, 1:8]
tfs_frame = frame[frame[:, 1] == 'TFS']
otf_frame = frame[frame[:, 1] == 'OTF']

tfs_num = 0
otf_num = 0
tfs_with_mass, tfs_size, tfs_total_mass, tfs_num = mass_handler(tfs_frame)
otf_with_mass, otf_size, otf_total_mass, otf_num = mass_handler(otf_frame)

tfs_com = mass_xyz(tfs_with_mass, tfs_total_mass, tfs_num, tfs_size)
otf_com = mass_xyz(otf_with_mass, otf_total_mass, otf_num, otf_size)

print tfs_com
print "\n------------------------\n"
print otf_com
 
