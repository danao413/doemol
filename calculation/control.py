from ase.io import write
import math
import numpy as np
import os
from pymatgen.io.cif import CifParser
from pymatgen.io.ase import AseAtomsAdaptor

path_to_species = '/projects/PTLearnPhoto/rtom/programs/FHI-aims/species_defaults'
species_default = 'light'
cif_dir = 'XXXI-2_structures_100'
batch_dir = 'mbd_light_batch_calc'

cif_dir = os.path.abspath(cif_dir)
batch_dir = os.path.abspath(batch_dir)
species_path = os.path.join(path_to_species, species_default)
cif_list = os.listdir(cif_dir)

def check_dir(dir_name):
	if not os.path.isdir(dir_name):
		os.mkdir(dir_name)

def unique_preserve_order(array):
	_, idx = np.unique(array, return_index=True)
	array = array[np.sort(idx)]
	return array

def write_control(struct, species_path, settings):
	#get the unique elements in the structure
	symbols = np.array(struct.get_chemical_symbols())
	symbols = unique_preserve_order(symbols)
	#get the corresponding atomic numbers
	ase_element = struct.get_atomic_numbers()
	element = []
	for ele in ase_element:
		if ele < 10:
			new_ele = str(0) + str(ele)
		else:
			new_ele = str(ele)
		element.append(new_ele)
	element = np.array(element)
	element = unique_preserve_order(element)
	control_file = 'control.in'
	with open(control_file, 'w') as f:
		for setting, value in settings.items():
			f.write('{}    {} \n'.format(setting, value))
	f.close()
	#create a list (species_list) of "species paths" from the elements and symbols in the array
	element = list(element)
	symbols = list(symbols)
	print(species_path)
	#TODO: fix this so that the elements correspond with the right symbols, the list is getting sorted somewhere... maybe in unique call
	species_list = []
	for ele, symbol in zip(element, symbols):
		species_list.append('{}_{}_default'.format(ele, symbol))
	species_list = [os.path.join(species_path, specie) for specie in species_list]
	for specie in species_list:
		with open(specie, 'r') as first_file, open(control_file, 'a') as second_file:
			for line in first_file:
				second_file.write(line)

def k_grid_25(struct):
	#can pass in the pymatgen lattice vectors and get a clear representation of the lattice vectors
	lattice = struct.lattice
	lv = list(lattice.abc)
	k_grid = [math.ceil(25/vector) for vector in lv]
	k_grid = '{} {} {}'.format(k_grid[0], k_grid[1], k_grid[2])
	return k_grid

#make a dictionary of the control settings for relaxed settings
relaxed_settings = {'xc': 'pbe',
					'spin': 'none',
					'relativistic': 'atomic_zora scalar',
					'charge': 0,
					'occupation_type': 'gaussian 0.01',
					'mixer': 'pulay',
					'n_max_pulay': 8,
					'charge_mix_param': 0.02, 
					'sc_accuracy_rho': 0.0001, 
					'sc_accuracy_eev': 0.01,
					'sc_accuracy_etot': 1e-06, 
					'sc_iter_limit': 10000, 
					'KS_method': 'parallel', 
					'empty_states': 6, 
					'basis_threshold': 1e-05,
					'k_grid': '3 3 3', 
					'sc_accuracy_forces': 0.0001, 
					'relax_geometry': 'trm 1e-2', 
					'relax_unit_cell': 'full',
					 'hessian_to_restart_geometry': '.false.', 
					 'harmonic_length_scale': 0.01, 
					 'energy_tolerance': 0.0005, 
					 'many_body_dispersion': ''}

#TODO: add the settings for this dictionary
SPE_settings = {'xc': 'pbe',
					'spin': 'none',
					'relativistic': 'atomic_zora scalar',
					'charge': 0,
					'occupation_type': 'gaussian 0.01',
					'mixer': 'pulay',
					'n_max_pulay': 8,
					'charge_mix_param': 0.02, 
					'sc_accuracy_rho': 0.0001, 
					'sc_accuracy_eev': 0.01,
					'sc_accuracy_etot': 1e-06, 
					'sc_iter_limit': 10000, 
					'KS_method': 'parallel', 
					'empty_states': 6, 
					'basis_threshold': 1e-05,
					'k_grid': '3 3 3', 
					'sc_accuracy_forces': 0.0001,  
					 'many_body_dispersion': ''}

check_dir(batch_dir)
os.chdir(cif_dir)

for cif in cif_list:
	cif_struct = CifParser(cif).get_structures()[0]
	ase_struct = AseAtomsAdaptor.get_atoms(cif_struct)
	struct_name = cif.split('.')[0]
	temp_path = os.path.join(batch_dir, struct_name)
	check_dir(temp_path)
	os.chdir(temp_path)
	write('geometry.in', ase_struct, format='aims')
	write('{}.in'.format(struct_name), ase_struct, format='aims')
	relaxed_settings['k_grid'] = k_grid_25(cif_struct)
	write_control(ase_struct, species_path, relaxed_settings)
	os.chdir(cif_dir)
