#generate mesh.conf file for phonopy
from ase.io import read
import numpy as np

#write the mesh.conf file
def write_mesh_conf(file_name, dim):
	#read in the supercell file and get the chemical symbols
	geo = read('geometry.in.supercell', format='aims')
	atom_name = geo.get_chemical_symbols()

	sc = geo.get_cell()
	lv = sc.lengths()

	mesh = 50 / lv + 1
	mesh = mesh.astype(int)

	with open(file_name, 'w') as f:
		f.write('DIM = {} {} {}\n'.format(dim, dim, dim))
		f.write('ATOM_NAME = {}\n'.format(atom_name))
		f.write('MP = {} {} {}\n'.format(mesh[0], mesh[1], mesh[2]))
		f.write('CUTOFF_FREQUENCY = -10\n')
		f.write('TPROP = .TRUE.\n')
		f.write('PRETEND_REAL = .TRUE.\n')
	f.close()
