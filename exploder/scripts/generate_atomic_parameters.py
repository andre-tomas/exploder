import os, sys
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.special as sc
import argparse

""" 
Code to parse atomic models from cretin
and create time-resolved data for Monte Carlo/molecular dynamics simulations.
Currenty, the code only works for atoms which are described by K, L, M shells (3).
Heavier atoms containing more shells (>3) can easily be incorporated, but requires
one to also change the GROMACS code to facilitate this. 
 
	Inputs: 
			photon_energy
	Outputs:	
			Rates for different processes. 
"""



def get_rate_data(cretin_atom_path):
	# template code for running rate equations with monte carlo
	fle = cretin_atom_path#sys.argv[1]

	#f = open("Kr_hyd.z01", "r")
	f = open(fle, "r")
	lines = f.readlines()
	atomic_state = []
	energy_level = []
	energy_level_dict = {}
	statistical_weight = []
	electronic_occupation = {}
	transition_dict = {}
	auger_transition_dict = {}
	phxs_transition_dict = {}
	phxs_transition_dict_inverse = {}
	fluorescence_transition_dict = {}
	statistical_weight_dictionary = {}
	collisional_ionization_dict = {}
	collisional_excitation_dict = {}

	atom_data = {}

	for line in lines:
		if "enot" in line:
			enot = float(line.split()[3])

		
		if "elev" in line:
			atomic_state.append([int(line.split()[len(line.split())-3]),int(line.split()[len(line.split())-2])])
			if len(line.split())==9:
				electronic_occupation[str([int(line.split()[1]),int(line.split()[2])])]=[int(line.split()[len(line.split())-3]),int(line.split()[len(line.split())-2]),0]
			else:
				electronic_occupation[str([int(line.split()[1]),int(line.split()[2])])]=[int(line.split()[len(line.split())-4]),int(line.split()[len(line.split())-3]),int(line.split()[len(line.split())-2])]
			statistical_weight.append(float(line.split()[4]))
			energy_level.append(float(line.split()[5]))
			if float(line.split()[5])==0:
				energy_level_dict[str([int(line.split()[1]),int(line.split()[2])])]=np.abs(enot) # remove np.abs later
			else:
				energy_level_dict[str([int(line.split()[1]),int(line.split()[2])])]=np.abs(enot-float(line.split()[5])) # remove np.abs later
				
			transition_dict[str([int(line.split()[1]),int(line.split()[2])])]=[]
			auger_transition_dict[str([int(line.split()[1]),int(line.split()[2])])]=[]
			phxs_transition_dict[str([int(line.split()[1]),int(line.split()[2])])]=[]
			phxs_transition_dict_inverse[str([int(line.split()[1]),int(line.split()[2])])]=[]
			fluorescence_transition_dict[str([int(line.split()[1]),int(line.split()[2])])]=[]
			statistical_weight_dictionary[str([int(line.split()[1]),int(line.split()[2])])]=float(line.split()[4])
			collisional_ionization_dict[str([int(line.split()[1]),int(line.split()[2])])]=[]
			collisional_excitation_dict[str([int(line.split()[1]),int(line.split()[2])])]=[]

	transition = []
	parameters = []

	for i, line in enumerate(lines):
		if "phis" in line:
			for j in range(i+1, 10000):
				if "end" in lines[j]:
					break
				if len(lines[j].split())>3:
					transition.append([[int(lines[j].split()[1]),int(lines[j].split()[2])], [int(lines[j].split()[3]),int(lines[j].split()[4])]])
					transition_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([int(lines[j].split()[3]),int(lines[j].split()[4])])
					
					parameters.append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])
					transition_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])


	auger_rates = [] # 1/sec
	auger_transition = []

	for i, line in enumerate(lines):
		if "augxs" in line:
			for j in range(i+1, 10000):
				if "end" in lines[j]:
					break
				if len(lines[j].split())>1:	
					lines[j].split()[1]
					auger_transition.append([[int(lines[j].split()[1]),int(lines[j].split()[2])],[int(lines[j].split()[3]),int(lines[j].split()[4])]])		
					auger_rates.append(float(lines[j].split()[len(lines[j].split())-1]))
					auger_transition_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([int(lines[j].split()[3]),int(lines[j].split()[4])])
					auger_transition_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append(float(lines[j].split()[len(lines[j].split())-1]))

	phxs_parameters = []
	for i, line in enumerate(lines):
		if "phxs" in line:
			for j in range(i+1, 10000):
				if "end" in lines[j]:
					break
				if len(lines[j].split())>3:
					phxs_transition_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([int(lines[j].split()[3]),int(lines[j].split()[4])])
					phxs_transition_dict_inverse[str([int(lines[j].split()[3]),int(lines[j].split()[4])])].append([int(lines[j].split()[1]),int(lines[j].split()[2])])
					phxs_parameters.append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])
					phxs_transition_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])
					phxs_transition_dict_inverse[str([int(lines[j].split()[3]),int(lines[j].split()[4])])].append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])


	iso1, i1 = 0, 0
	for i, line in enumerate(lines):
		if "augis" in line:
			for j in range(i+1, 10000):
				if "end" in lines[j]:
					break
				if len(lines[j].split())>0:
					if "d"== lines[j].split()[0]: #and not "done" in lines[j] and not "rad" in lines[j]:
						iso1, i1 = int(lines[j].split()[1]), int(lines[j].split()[4])
						
				if "rad" in lines[j]:
					fluorescence_transition_dict[str([iso1, i1])].append([int(lines[j].split()[1]),int(lines[j].split()[2])])
					fluorescence_transition_dict[str([iso1, i1])].append([float(x) for x in lines[j].split()[3:len(lines[j].split())]])

	
	collisional_excitation_parameters = []
	for i, line in enumerate(lines):
		if "colex2" in line:
			for j in range(i+1, 10000):
				if "end" in lines[j]:
					break
				if len(lines[j].split())>3:
					collisional_excitation_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([int(lines[j].split()[3]),int(lines[j].split()[4])])
					collisional_excitation_parameters.append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])
					collisional_excitation_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])



	collisional_parameters = []
	for i, line in enumerate(lines):
		if "sampson" in line:
			for j in range(i+1, 10000):
				if "end" in lines[j]:
					break
				if len(lines[j].split())>3:
					collisional_ionization_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([int(lines[j].split()[3]),int(lines[j].split()[4])])

					collisional_parameters.append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])
					collisional_ionization_dict[str([int(lines[j].split()[1]),int(lines[j].split()[2])])].append([float(x) for x in lines[j].split()[5:len(lines[j].split())]])

	return atomic_state, electronic_occupation, statistical_weight, energy_level, transition_dict, transition_dict, auger_transition_dict, phxs_transition_dict, phxs_transition_dict_inverse, fluorescence_transition_dict, statistical_weight_dictionary, energy_level_dict, collisional_ionization_dict, collisional_excitation_dict


def get_molecule_data(xyz_file):
	f_molecule = open(xyz_file, "r")
	lines = f_molecule.readlines()

	atomic_species = []
	atoms = []
	xyz = []
	number_of_atoms=0
	for line in lines:
		add = True
		try:
			if len(line.split())==1:
				number_of_atoms = int(line.split()[0])
			else:
				atoms.append(line.split()[0])
				xyz.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
				for element in atomic_species:
					if element == line.split()[0]:
						add = False
				if add == True:
					atomic_species.append(line.split()[0])
		except:
			continue

	if 'H' not in atomic_species:
		atomic_species.append('H')

	return atoms, xyz,number_of_atoms,atomic_species

def gaussian(v, mu, sigma):
    """v-freq, mu-center, sigma-width"""
    return (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((v-mu)/sigma)**2)

def lorentzian(v, mu, gamma):
    """v-freq, mu-center, gamma-1/height"""
    return (1/np.pi)*(gamma/((v-mu)**2 + gamma**2))

def photoexcitation_cross_section(v, f, photon_frequency):
	sigma = 1
	c = 29979245800 # [cm / s]
	e = 1.60217663*1e-19 
	m = 9.1093837*1e-31
	return ((np.pi*e**2*f)/(m*c))*gaussian(v, photon_frequency, sigma)

def get_photon_frequency(energy):
	energy = energy/(6.24150907*1e18)
	h = 6.62607015*1e-34 # Plankc's constant
	return energy/h

def get_photon_photon_energy(frequency):
	h = 6.62607015*1e-34 # Plankc's constant
	return h*frequency*6.24150907*1e18

def calculate_photoexcitation_rate(state, fluence, photon_energy):

	final_states = []
	parameters = []
	photoexcitation_rates = []
	indx = 0
	c = 29979245800 # [cm / s]
	e = 1.60217663*1e-19 
	m = 9.1093837*1e-31
	for i, final_state in enumerate(phxs_transition_dict[state]):
	
		if i%2==0:
			final_states.append(final_state)
		else:
			if len(phxs_transition_dict[state][i]) == 3:
				f, LAMBDA, width = phxs_transition_dict[state][i][0], phxs_transition_dict[state][i][1], phxs_transition_dict[state][i][2]
			if len(phxs_transition_dict[state][i]) == 2:
				f, LAMBDA, width = phxs_transition_dict[state][i][0], phxs_transition_dict[state][i][1], 0
			LAMBDA = LAMBDA*1e-10
			photon_frequency = c/LAMBDA
			v = get_photon_frequency(photon_energy)
			v = photon_frequency
			phxs_cs = photoexcitation_cross_section(v, f, photon_frequency)
			photoexcitation_rates.append(fluence*phxs_cs)
	

	return final_states, photoexcitation_rates

def calculate_fluorescence_rate(state):
	""" Use same oscillator strength as for photoexcitation. """

	c = 29979245800 # [cm / s]
	e = 1.60217663*1e-19 
	m = 9.1093837*1e-31
	final_states = []
	parameters = []
	fluorescence_rates = []
	indx = 0
	for i, final_state in enumerate(phxs_transition_dict_inverse[state]):
		if i%2==0:
			final_states.append(final_state)
		else:
			if len(phxs_transition_dict_inverse[state][i]) == 3:
				f, LAMBDA, width = phxs_transition_dict_inverse[state][i][0], phxs_transition_dict_inverse[state][i][1], phxs_transition_dict_inverse[state][i][2]
			if len(phxs_transition_dict_inverse[state][i]) == 2:
				f, LAMBDA, width = phxs_transition_dict_inverse[state][i][0], phxs_transition_dict_inverse[state][i][1], 0
			LAMBDA = LAMBDA*1e-10
			v = c/LAMBDA
			gi = statistical_weight_dictionary[state]
			gf = statistical_weight_dictionary[str(final_states[indx])]
			indx+=1
			fluorescence_rate = (8*np.pi**2*e**2*v**2/(m*c**3))*(gi/gf)*f
			fluorescence_rates.append(fluorescence_rate)
			

	return final_states, fluorescence_rates


def photoionization_crossection(a0, a1, a2, a3, n, de, emin, emax, energy):
	#a0, a1, a2, a3, n, de, emin, emax  = parameters[5]
	#energy = np.linspace(emin, emax, 10000)

	if emin<= np.min(energy) and np.max(energy)<=emax:

		#energy =8000
		b = np.divide(energy,de)
		pblog = a0 + a1*np.log(b) + a2*np.log(b)**2 + a3*np.log(b)**3  
		pb = np.exp(pblog)
		sigma = 1e-18*n*pb*(13.606/de)
		#sigma = 3*1e28*1e-6*1e-2*1e-2*1e-18*n*pb*(13.606/de) # unit when comparing to other references
		return sigma

	else:
		#Energy out of bonds. Setting photoionization crossection to zero.
		return 0

		
def get_photoionization_cs(energy, state):

	phis_prob = []

	for j in range(len(transition_dict[state])):
		if (j % 2) == 0:
			continue
		else:
			parameters = transition_dict[state][j]
			a0, a1, a2, a3, n, de, emin, emax  = parameters
			phis_prob.append(photoionization_crossection(a0, a1, a2, a3, n, de, emin, emax, energy))

	return phis_prob

def get_electronic_state(state):
	return electronic_occupation[state]

def get_transition_electronic_states(state):
	electronic_states = []
	for j in range(len(transition_dict[state])):
			if (j % 2) == 0:
				estate = transition_dict[state][j]
				electronic_states.append(electronic_occupation[str(estate)])
			else:
				continue
	return electronic_states

def get_transition_electronic_states_auger(state):
	electronic_states = []
	for j in range(len(auger_transition_dict[state])):
			if (j % 2) == 0:
				estate = auger_transition_dict[state][j]
				electronic_states.append(electronic_occupation[str(estate)])
			else:
				continue
	return electronic_states

def get_transition_electronic_states_phxs(states):
	electronic_states = []
	for j in range(len(states)):
			estate = states[j]
			electronic_states.append(electronic_occupation[str(estate)])
		
	return electronic_states

def get_name_of_electron_occupation(K,L,M,electronic_occupation):
	for estate in electronic_occupation:
		if K == (int(electronic_occupation[estate][0])) and L == (int(electronic_occupation[estate][1])) and M == (int(electronic_occupation[estate][2])):
			return estate

def get_auger_rates(state):
	auger_rates_t = []
	for j in range(len(auger_transition_dict[(state)])):
		if (j % 2) == 0:
			continue
		else:
			values = auger_transition_dict[(state)][j]
			auger_rates_t.append(values)

	return auger_rates_t

def photoionization_probability(crossection, fluence, dt):
	return np.multiply(crossection, fluence*dt)

def photoionization_rates(crossection, fluence):
	return np.multiply(crossection, fluence)

def concatenate_transitions(photoionization, auger):
	
	if len(photoionization) == 0 and len(auger) == 0:
		return []
	elif len(photoionization) == 0 and len(auger)>0:
		return auger
	elif len(photoionization) > 0 and len(auger)==0:
		return photoionization
	else:
		return np.concatenate((photoionization, auger))

def run_monte_carlo_rate_step(rate_list, uQ):
	for j in range(1, len(rate_list)):
		if rate_list[j-1] < uQ <= rate_list[j]:
			return rate_list[j]
	return 0

def randomise(N):
	random_numbers = []
	for i in range(N):
		u = random.uniform(0,1)
		random_numbers.append(u)
	return random_numbers

def collisonal_crossection(data, tev):
	c0, c1, c2, c3, s0 = data[0], data[1], data[2], data[3], data[4]
	delta_E = 5
	d = delta_E/tev
	#tev = 5
	Ei = sc.expi(b)
	q = s0*(c0 + c1*np.sqrt(np.pi*b)*sc.erfc(b)*np.exp(b**2) + (c2*b + c3)*Ei*np.exp(b))
	sigma = (1.09*1e-6*q*np.exp(-b))/(delta_E*np.sqrt(tev))
	return sigma

def get_collisonal_rates(state, tev):
	collisional_rates_t = []

	for j in range(len(collisional_ionization_dict[(state)])):
		if (j % 2) == 0:
			continue
		else:
			values = collisional_ionization_dict[(state)][j]
			rate = collisonal_crossection(values, tev)
			collisional_rates_t.append(rate)

	return collisional_rates_t

def get_transition_electronic_states_collisional(state):
	electronic_states = []
	for j in range(len(collisional_ionization_dict[state])):
			if (j % 2) == 0:
				estate = collisional_ionization_dict[state][j]
				electronic_states.append(electronic_occupation[str(estate)])
			else:
				continue
	return electronic_states
	
def monte_carlo_electronic_transition_XMDYN(state,photon_energy):

	phis_cs = get_photoionization_cs(photon_energy,state)
	nonzero_inds=np.where(np.array(phis_cs)>0)[0]
	phis_cs=np.array(phis_cs)[nonzero_inds]	
	phis_rate = photoionization_rates(phis_cs, fluence)
	phis_estates = get_transition_electronic_states(state)
	phis_estates=np.array(phis_estates)[nonzero_inds]
	
	auger_rates_t = get_auger_rates(state)
	nonzero_inds=np.where(np.array(auger_rates_t)>0)[0]
	auger_rates_t=(np.array(auger_rates_t)[nonzero_inds])
	auger_estates = get_transition_electronic_states_auger(state)
	auger_estates=np.array(auger_estates)[nonzero_inds]

	collisional = False
	if collisional == True:
		collisonal_rates_t = get_collisonal_rates(state, tev)
		nonzero_inds=np.where(np.array(collisonal_rates_t)>0)[0]
		collisonal_rates_t = (np.array(collisonal_rates_t)[nonzero_inds])
		collisional_estates = get_transition_electronic_states_collisional(state)
		collisional_estates = np.array(collisional_estates)[nonzero_inds]

	final_states, fluorescence_rates=calculate_fluorescence_rate(state)
	nonzero_inds=np.where(np.array(fluorescence_rates)>0)[0]
	fluorescence_rates=np.array(fluorescence_rates)[nonzero_inds]	
	fluorescence_estates = get_transition_electronic_states_phxs(final_states)
	
	#photon_energy = 653.7113303908512
	final_states, phxs_rates = calculate_photoexcitation_rate(state, 1, photon_energy)
	nonzero_inds=np.where(np.array(phxs_rates)>0)[0]
	phxs_rates=np.array(phxs_rates)[nonzero_inds]	
	phxs_estates = get_transition_electronic_states_phxs(final_states)

	transition_list = concatenate_transitions(phis_estates, auger_estates)
	transition_list= concatenate_transitions(transition_list, fluorescence_estates)

	if len(transition_list)==0:
		return [], []
	# Rate calculations
	total_rate = np.sum(phis_rate) + np.sum(auger_rates_t) + np.sum(phxs_rates) + np.sum(fluorescence_rates)# + np.sum(collisonal_rates_t)
	
	# 1st of Monte Carlo algorithm (finding transition state)
	rate_list = np.concatenate((phis_rate,auger_rates_t))
	rate_list = np.concatenate((rate_list,fluorescence_rates))
	rate_list = np.array(rate_list)
	transition_list = np.array(transition_list)
	inds = rate_list.argsort()
	rate_list = rate_list[inds]
	sorted_transition = transition_list[inds]

	return rate_list, sorted_transition

def initialize_atomic_states(atoms):
	atom_to_state = {}
	atom_paths = {}
	atom_to_state["H"] = [1,0,0]
	atom_to_state["O"] = [2,6,0]
	atom_to_state["N"] = [2,5,0]
	atom_to_state["C"] = [2,4,0]
	atom_to_state["P"] = [2,8,5]
	atom_to_state["S"] = [2,8,6]
	atom_to_state["Fe"] = [2,8,16]
	atomic_numbers = []
	Z = {}
	Z["O"] = 8
	Z["H"] = 1
	Z["S"] = 16
	Z["P"] = 15
	Z["N"] = 7
	Z["C"] = 6
	Z["Fe"] = 26
	atom_paths["O"]  = f"{main_path}/Atomic_model/oxygen/O_SCH_model.z01"
	atom_paths["H"]  = f"{main_path}/Atomic_model/hydrogen/H_SCH_model.z01"
	atom_paths["C"]  = f"{main_path}/Atomic_model/carbon/C_SCH_model.z01"
	atom_paths["N"]  = f"{main_path}/Atomic_model/nitrogen/N_SCH_model.z01"
	atom_paths["P"]  = f"{main_path}/Atomic_model/phosphor/P_SCH_model.z01"
	atom_paths["S"]  = f"{main_path}/Atomic_model/sulfur/S_SCH_model.z01"
	atom_paths["Fe"] = f"{main_path}/Atomic_model/iron/Fe_SCH_model.z01"

	atomic_states = []
	for atom in atoms:
		atomic_state, electronic_occupation, statistical_weight, energy_level, transition_dict, transition_dict, auger_transition_dict, phxs_transition_dict, phxs_transition_dict_inverse, fluorescence_transition_dict, statistical_weight_dictionary, energy_level_dict, collisional_ionization_dict, collisional_excitation_dict = get_rate_data(atom_paths[atom])

		K, L, M= atom_to_state[atom][0], atom_to_state[atom][1], atom_to_state[atom][2]
		atomic_states.append(get_name_of_electron_occupation(K,L,M,electronic_occupation))
		atomic_numbers.append(Z[atom])

	return atomic_states, atomic_numbers, atom_paths

def intialize_cretin_data():
	atoms, xyz, number_of_atoms, atomic_species = get_molecule_data(xyz_file)

	CRETIN_DATA = {}
	for j, atom in enumerate(atomic_species):
		cretin_atom_path = atom_paths[atom]
		atomic_state, electronic_occupation, statistical_weight, energy_level, transition_dict, transition_dict, auger_transition_dict, phxs_transition_dict, phxs_transition_dict_inverse, fluorescence_transition_dict, statistical_weight_dictionary, energy_level_dict, collisional_ionization_dict, collisional_excitation_dict = get_rate_data(cretin_atom_path)
		CRETIN_DATA[atom] = get_rate_data(cretin_atom_path)
	return CRETIN_DATA

def get_atom_species():
	unique_atom_list = []

	for atom in atoms:
		add=True
		for atm2 in unique_atom_list:
			if atom == atm2:
				add = False
				break
		if add ==True:
			unique_atom_list.append(atom)
	return unique_atom_list

def convert_pdb_to_xyz(pdb_file, output_name):
	""" Function to convert pdb file to xyz using obabel. """
	os.system(f"obabel -i pdb {src}/{pdb_file} -o xyz -O {dst}/{output_name}.xyz")
	return f"{dst}/{output_name}.xyz"

### Modify output files ###

# Separation = 4 for rates and 8 for colls
def add_delimiter_to_file(filename,separation):
	flag = ("_O" in filename) and ("coll" in filename)

	with open(filename, 'r') as file:
		lines = file.readlines()

	modified_lines = []

	for line in lines:
		numbers = line.split()
		if flag:
			print(line)
   	 
    	# Add the first three numbers without a semicolon
		modified_line = ' '.join(numbers[:3])
    
    	# Add the remaining numbers in groups of four with semicolons
		for i in range(3, len(numbers), separation):
			modified_line += ' ; ' + ' '.join(numbers[i:i+separation])

		modified_lines.append(modified_line + ' ;\n')

	with open(filename,'w') as file:
		for line in modified_lines:

			file.write(line)


def extract_unique_elements_from_pdb(file_path):
    """
    Extracts unique elements from a PDB file.
    
    Args:
    file_path (str): Path to the PDB file

    Returns:
    list: A list of unique elements in the PDB file
    """
    unique_elements = set()

    # Read the PDB file
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract the element symbol, which is typically found at these positions
                # The element symbol is right-aligned within columns 77-78
                element = line[76:78].strip()
                unique_elements.add(element)

    return list(unique_elements)





	

### END OF FUNCTIONS ###
parser = argparse.ArgumentParser()
parser.add_argument("photon_energy", help="Photon energy in eV.")
parser.add_argument("sample_file", help="Name of input (PDB or XYZ).")
parser.add_argument("src", help="PDB/XYZ source path.")
parser.add_argument("dst", help="Output destination path.")
args = parser.parse_args()

photon_energy = args.photon_energy
sample_file = args.sample_file
src = args.src
dst = args.dst

run_path = os.getcwd()
main_path  = os.path.realpath(os.path.dirname(__file__))

# Get the absolute paths
os.chdir(src)
src = os.getcwd()
os.chdir(run_path)
os.chdir(dst)
dst = os.getcwd()


try:
	os.mkdir("Atomic_data")
	os.chdir("Atomic_data")
except:
	os.system("rm -rf Atomic_data")
	os.mkdir("Atomic_data")
	os.chdir("Atomic_data")


if "pdb" in sample_file:
	pdb_file = sample_file
	xyz_file = convert_pdb_to_xyz(pdb_file, "sample")
else:
	xyz_file = f"{src}/{sample_file}" 

atoms, xyz, number_of_atoms, atomic_species = get_molecule_data(xyz_file) # Get info regarding atom species in the system and coordinates


atomic_states, atomic_numbers, atom_paths = initialize_atomic_states(atoms) # initializing electronic states and possible transitions
CRETIN_DATA = intialize_cretin_data()

fluence = 1 # Not used anymore, we use the fluence in GROMACS
photon_energy = float(args.photon_energy) # 600, Photon energy, this is the only parameter we change. Will change the cross section for photoprocesses
atom = atoms[1]
state = atomic_states[1]

unique_atom_list = get_atom_species()
if 'H' not in unique_atom_list: # We always wany H but it is sometimes omitted from pdb/xyz
		unique_atom_list.append('H')




atomic_state, electronic_occupation, statistical_weight, energy_level, transition_dict, transition_dict, auger_transition_dict, phxs_transition_dict, phxs_transition_dict_inverse, fluorescence_transition_dict, statistical_weight_dictionary, energy_level_dict, collisional_ionization_dict, collisional_excitation_dict= CRETIN_DATA[atom]

recombination = True 
if recombination == False:
	# No electronic recombination.

	for atom in unique_atom_list:



		atomic_state, electronic_occupation, statistical_weight, energy_level, transition_dict, transition_dict, auger_transition_dict, phxs_transition_dict, phxs_transition_dict_inverse, fluorescence_transition_dict, statistical_weight_dictionary, energy_level_dict, collisional_ionization_dict, collisional_excitation_dict = CRETIN_DATA[atom]
		fcoll = open("collisional_parameters_{}.txt".format(atom.upper()), "w")

		for state in electronic_occupation:
			estate = get_electronic_state(state)
			fcoll.writelines("{0} {1} {2} ".format(estate[0], estate[1], estate[2]))

			for j in range(len(collisional_ionization_dict[(state)])):
				if (j % 2) == 0:
					estate = collisional_ionization_dict[(state)][j]
					estate = (get_electronic_state(str(estate)))
				else:
					values = collisional_ionization_dict[(state)][j]
					fcoll.writelines("{0} {1} {2} {3} {4} {5} {6} {7} ".format(estate[0], estate[1], estate[2],values[0],values[1],values[2],values[3],values[4]))
			fcoll.writelines("\n")
		fcoll.close()

else:

	for atom in unique_atom_list:
		atomic_state, electronic_occupation, statistical_weight, energy_level, transition_dict, transition_dict, auger_transition_dict, phxs_transition_dict, phxs_transition_dict_inverse, fluorescence_transition_dict, statistical_weight_dictionary, energy_level_dict, collisional_ionization_dict, collisional_excitation_dict = CRETIN_DATA[atom]
		fcoll = open("collisional_parameters_{}.txt".format(atom.upper()), "w")

		for state in electronic_occupation:
			estate = get_electronic_state(state)

			fcoll.writelines("{0} {1} {2} ; ".format(estate[0], estate[1], estate[2]))
			# collisional ionization
			for j in range(len(collisional_ionization_dict[(state)])):
				if (j % 2) == 0:
					estate = collisional_ionization_dict[(state)][j]
					estate = (get_electronic_state(str(estate)))
				else:
					values = collisional_ionization_dict[(state)][j]
					fcoll.writelines("{0} {1} {2} {3} {4} {5} {6} {7} ; ".format(estate[0], estate[1], estate[2],values[0],values[1],values[2],values[3],values[4]))
			for state2 in electronic_occupation:
				if state!=state2:
				
					for j in range(len(collisional_ionization_dict[(state2)])):
						if (j % 2) == 0:
							

							ESTATE = collisional_ionization_dict[(state2)][j]
							if str(ESTATE) == str(state):
								values = collisional_ionization_dict[(state2)][j+1]
								recomb_state = get_electronic_state(state2)	
								fcoll.writelines("{0} {1} {2} {3} {4} {5} {6} {7} ; ".format(recomb_state[0], recomb_state[1], recomb_state[2],values[0],values[1],values[2],values[3],values[4]))
						
			for j in range(len(collisional_excitation_dict[(state)])):
				if (j % 2) == 0:
					estate = collisional_excitation_dict[(state)][j]
					estate = (get_electronic_state(str(estate)))
				else:
					values = collisional_excitation_dict[(state)][j]
					fcoll.writelines("{0} {1} {2} {3} {4} {5} {6} {7} ; ".format(estate[0], estate[1], estate[2],values[0],values[1],values[2],values[3],values[4]))

			fcoll.writelines("\n")

for atom in unique_atom_list:

	f2 = open("rate_transitions_{}.txt".format(atom.upper()), "w")
	
	atomic_state, electronic_occupation, statistical_weight, energy_level, transition_dict, transition_dict, auger_transition_dict, phxs_transition_dict, phxs_transition_dict_inverse, fluorescence_transition_dict, statistical_weight_dictionary, energy_level_dict, collisional_ionization_dict, collisional_excitation_dict = CRETIN_DATA[atom]

	for state in electronic_occupation:

		initial_state = (get_electronic_state(state))
		if initial_state[0] == 0 and initial_state[1] == 0 and initial_state[2]==0:
			continue

		rate_list, sorted_transition = monte_carlo_electronic_transition_XMDYN(state,photon_energy)

		if (len(rate_list)==0):
			f2.writelines("{0} {1} {2} ; 0 0 0 0.000 ; \n".format(initial_state[0], initial_state[1], initial_state[2]))
			continue

		stringval = "{0} {1} {2} ;".format(initial_state[0], initial_state[1], initial_state[2])
		for transition, rate in zip(sorted_transition, rate_list):
			stringval += "{0} {1} {2} {3} ;".format(transition[0], transition[1], transition[2], rate) 

		f2.writelines(stringval + "\n")

	f2.writelines("0 0 0 ; 0 0 1  0.000 ; " + "\n")
	f2.close()

	f2 = open("energy_levels_{}.txt".format(atom.upper()), "w")
	for state in electronic_occupation:

		estate = get_electronic_state(state)
		f2.writelines("{0} {1} {2} {3}\n".format(estate[0], estate[1], estate[2], energy_level_dict[state]))

	f2.close()
	f2 = open("statistical_weight_{}.txt".format(atom.upper()), "w")
	for state in electronic_occupation:
		estate = get_electronic_state(state)
		f2.writelines("{0} {1} {2} {3}\n".format(estate[0], estate[1], estate[2], statistical_weight_dictionary[state]))
	f2.close()

## Modify files
'''
os.chdir(f"{dst}/Atomic_data")
files = os.listdir(os.getcwd())

for file in (file for file in files if "rate" in file):
	add_delimiter_to_file(file,4)
for file in (file for file in files if "coll" in file):
	print(file)
	add_delimiter_to_file(file,8)
'''

#os.system(f"rm {dst}/sample.xyz")




