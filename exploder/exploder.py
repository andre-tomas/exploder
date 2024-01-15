import os,sys,shutil
import numpy as np
import shlex 
from subprocess import call
import subprocess as sp
import MDAnalysis as md
import matplotlib.pyplot as plt
import time
import pickle
import h5py
import datetime

# MC EDITION

'''
  ________ __    ___      ___ ____  ____ ___           __ ___________ __     ______   _____  ___   
 /"       )" \  |"  \    /"  ("  _||_ " |"  |         /""("     _   ")" \   /    " \ (\"   \|"  \  
(:   \___/||  |  \   \  //   |   (  ) : ||  |        /    )__/  \\__/||  | // ____  \|.\\   \    | 
 \___  \  |:  |  /\\  \/.    (:  |  | . ):  |       /' /\  \ \\_ /   |:  |/  /    ) :): \.   \\  | 
  __/  \\ |.  | |: \.        |\\ \__/ // \  |___   //  __'  \|.  |   |.  (: (____/ //|.  \    \. | 
 /" \   :)/\  |\|.  \    /:  |/\\ __ //\( \_|:  \ /   /  \\  \:  |   /\  |\        / |    \    \ | 
(_______/(__\_|_)___|\__/|___(__________)\_______|___/    \___)__|  (__\_|_)"_____/   \___|\____\
Simulation (name pending) is a class whose job is to integrate with GROMACS.
Each Simulation is connected to a .mdp (md parameter file) which it uses to run the sim.

'''

class Simulation():
	def __init__(self,
		path_to_mdp,
		gmx_path,
		forcefield,
		threads,
		ts,
		nsteps,
		gaussian_peak,
		num_photons,
		sigma,
		focal_diameter,
		photon_energy,
		do_ff=1,
		do_ct=1,
		do_debye=1,
		do_log=0,
		read_states=0,
		double=False,
		MD_log=1000):

		# mdp parameters
		self.ts = ts
		self.nsteps = nsteps
		self.gaussian_peak = gaussian_peak
		self.num_photons = num_photons
		self.sigma = sigma
		self.focal_diameter = focal_diameter
		self.photon_energy = photon_energy
		self.MD_log = MD_log


		# What do?
		self.do_ff = do_ff
		self.do_ct = do_ct
		self.do_debye=do_debye
		self.do_log=do_log
		self.read_states=read_states


		# other parameters
		self.forcefield = forcefield
		self.threads = threads

		# sim files
		self.mdp = path_to_mdp

		# Gromacs paths
		self.path 		= gmx_path
		self.grompp 	= f"{gmx_path}/grompp"
		self.mdrun 	= f"{gmx_path}/mdrun"
		self.pdb2gmx 	= f"{gmx_path}/pdb2gmx"
		self.editconf 	= f"{gmx_path}/edifconf"
		self.trjconv 	= f"{gmx_path}/trjconv"

		if double:
			self.grompp 	= f"{self.grompp}_d"
			self.mdrun 	= f"{self.mdrun}_d"
			self.pdb2gmx 	= f"{self.pdb2gmx}_d"
			self.editconf 	= f"{self.editconf}_d"
			self.trjconv 	= f"{self.trjconv}_d"
			


		self.configure()

	def __str__(self):
		return f"Simulation using:{self.path} \n with mdp: {self.mdp} \n runs {self.nsteps} with ts {self.ts*1000}fs"

	def __repr__(self):
		return self.__str__()

	def configure(self):
		with open(self.mdp,'r') as f:
			lines = f.readlines()

		with open(self.mdp,"w") as f:
			for line in lines:
				split_line =line.split()

				if "dt" in split_line:
					line = f"dt                       = {self.ts}; \n"
				if "nsteps" in split_line:
					line = f"nsteps                   = {self.nsteps}; \n"
				if "userreal1" in split_line: 
					line = f"userreal1				  = {self.gaussian_peak};\n"
				if "userreal2" in split_line: 
					line = f"userreal2				  = {self.num_photons}; \n"
				if "userreal3" in split_line:
					line = f"userreal3				  = {self.sigma}; \n"
				if "userreal4" in split_line:
					line = f"userreal4				  = {self.focal_diameter}; \n"
				if "userreal5" in split_line:
					line = f"userreal5				  = {self.photon_energy}; \n"



				if "userint1" in split_line:
					line = f"userint1				  = {self.do_ff}; \n"
				if "userint2" in split_line:
					line = f"userint2				  = {self.do_ct}; \n"
				if "userint3" in split_line:
					line = f"userint3				  = {self.do_debye}; \n"
				if "userint4" in split_line:
					line = f"userint4				  = {self.read_states}; \n"
				if "userint5" in split_line:
					line = f"userint5				  = {self.do_log}; \n"

				# Log stuff
				if "nstxout" in split_line: 
					line = f"nstxout					= {self.MD_log}; \n"
				if "nstvout" in split_line: 
					line = f"nstvout					= {self.MD_log}; \n"
				if "nstfout" in split_line: 
					line = f"nstfout					= {self.MD_log}; \n"
				if "nstlog" in split_line: 
					line = f"nstlog					= {self.MD_log}; \n"
				if "nstenergy" in split_line: 
					line = f"nstenergy					= {self.MD_log}; \n"
				if "xtc-precision" in split_line: 
					line = f"xtc-precision					= {self.MD_log}; \n"





				f.write(line)




	def gmx_pdb2gmx(self,input_file,output_file,rotate="0 0 0"):

		file = input_file.split('.')[0]
		ext  = input_file.split('.')[-1]


		cmd = f"{self.pdb2gmx} -f {input_file} -o pdb2gmx_intermediate.gro -water spce -ff {self.forcefield} -ignh"
		sp.Popen(shlex.split(cmd)).wait()

		cmd = f"{gmx_path}/editconf -f pdb2gmx_intermediate.gro -o {output_file} -rotate {rotate} -c -d 12 -bt cubic"
		call(shlex.split(cmd))

	def gmx_grompp(self,input_file, output_file,topol_file="topol.top"):

		cmd = f"{self.grompp} -f {self.mdp} -c {input_file} -p {topol_file} -o {output_file} -maxwarn 5"
		sp.Popen(shlex.split(cmd)).wait()

# The parallel keyword "par" does not paralallize each simulation! but rather lets us run multiple simulationS at each time.

	def gmx_mdrun(self,tpr_name,ionize=False,par=False):

		if tpr_name.split(".")[-1] == "tpr":
			tpr_name = tpr_name.split(".")[0]

		if par:
			cmd = f"{self.mdrun} -deffnm {tpr_name} -v -nt 1 -rdd 1000"
		else:
			cmd = f"{self.mdrun} -deffnm {tpr_name} -v -nt {self.threads} -rdd 1000"


		if ionize:
			cmd = f"{cmd} -ionize"

		if par: # Must be last
			cmd = f"{cmd}"
			return sp.Popen(shlex.split(cmd))
		
		else: 
			call(shlex.split(cmd))

	def gmx_grunpp(self,grompp_input,name,ionize=False,par=False,topol_file ="topol.top"):
		self.gmx_grompp(grompp_input,name,topol_file)
		if par:
			p = self.gmx_mdrun(name,ionize,par)
			return p
		else:
			self.gmx_mdrun(name,ionize,par)

	def gmx_trjconv(self,trr,tpr,output_path,time,choice=0):
		cmd = f"echo '{choice}' | {self.trjconv} -f {trr} -s {tpr} -dump {time} -o {output_path}"  
		os.system(cmd)



"""
________________________________________________________________________________________________________
 ___      ___    ______   ________   _______ ___       
|"  \    /"  |  /    " \ |"      "\ /"     "|"  |      
 \   \  //   | // ____  \(.  ___  :|: ______)|  |      
 /\\  \/.    |/  /    ) :): \   ) ||\/    | |:  |      
|: \.        (: (____/ //(| (___\ ||// ___)_ \  |___   
|.  \    /:  |\        / |:       :|:      "( \_|:  \  
|___|\__/|___| \"_____/  (________/ \_______)\_______)


"""
class Model():

	def __init__(self,name,input_path,result_path,conf_sim,exp_sim1,exp_sim2):

		self.name 				= name
		self.input 				= input_path
		self.exp_sim1 	  = exp_sim1
		self.exp_sim2 		= exp_sim2
		self.conf_sim 	 	= conf_sim
		self.result_path 	= f"{result_path}/{self.name}"
		self.output_path 	= result_path

		#MC 
		self.atomic_data = f"{main_path}/Atomic_data"
		self.cretin_data = f"{main_path}/simulation_output"


		os.mkdir(self.result_path)
		os.chdir(self.result_path)

		self.simulation_path = f"{self.result_path}/simulation"

		os.mkdir(self.simulation_path)

	def __str__(self): 
		return f"Model: {self.name} \n"
	
	def __repr__(self):
		return self.__str__()


	def configuration_simulation(self):
		os.chdir(self.simulation_path)

		gro_file = f"{self.name}.gro"
		deffnm = f"{self.name}_conf"
		self.conf_sim.gmx_pdb2gmx(self.input,gro_file)
		self.conf_sim.gmx_grunpp(gro_file,deffnm)

		t_end = self.conf_sim.nsteps*self.conf_sim.ts
		t_start = np.round(t_end/3)




		for l,t in enumerate(np.linspace(t_start,t_end,extracted_frames)):	
			dump_output = f"{self.simulation_path}/{self.name}_exp_{l:05d}.gro"
			self.conf_sim.gmx_trjconv(f"{deffnm}.trr",f"{deffnm}.tpr",dump_output,t)

		self.exp_names = f"{self.name}_exp"

	'''
	Running in parallel does not mean each simulation is running with multiple cores, 
	but rather simultainous simulations are running. 
	
	'''
	def explosion_simulation(self):
		sim_dir = os.listdir(self.simulation_path)
		self.explode_files = sorted([x for x in sim_dir if self.exp_names in x])

		

		sweeps = extracted_frames/params.threads_exp
		if not sweeps.is_integer():
			print("Number of frames must be divisible by number of threads for parallell simulations.")
			print(f"Setting sweeps to {(extracted_frames//params.threads_exp)}")
			sweeps = extracted_frames//params.threads_exp
		sweeps = int(sweeps)



		for sweep in range(sweeps):
			print(f"Entering sweep {sweep+1}/{sweeps}")

			proc = [] 
			## Initial explosions steps use very small timestep
			for configuration in self.explode_files[params.threads_exp*sweep:(sweep+1)*params.threads_exp]:
				config_name = configuration.split('.')[0]
				os.chdir(self.simulation_path)
				os.mkdir(config_name)
				shutil.copyfile(configuration,f"{config_name}/{configuration}")
				os.chdir(config_name)
				os.mkdir("simulation_output") # Used to write electron dynamics logs

				generate_atomic_parameters(600,f"{self.name}.pdb",input_path,"./")
				
				with open(log_ext_path,'a') as f:
					log = f"Starting run {config_name} from {self.name} at time: {datetime.datetime.now()}.\n"

				p = self.exp_sim1.gmx_grunpp(configuration,config_name,ionize=True,par=True,topol_file="../topol.top")


				proc.append(p)

			for p in proc: p.wait()
			counts = np.ones(shape=len(proc))

			doneflag = [False for _ in range(len(proc))]

			# We keep running the sim until a threshold energy is reached.
			while not all(doneflag):
				proc = []
				for k,configuration in enumerate(self.explode_files[params.threads_exp*sweep:(sweep+1)*params.threads_exp]):
					config_name = configuration.split('.')[0]
					os.chdir(self.simulation_path)
					os.chdir(config_name)

					cmd = f"echo -ne 'Potential \n Total-Energy' | {self.exp_sim2.path}/g_energy -f {config_name}.edr"  
					os.system(cmd)		
					t,pot,tot =  read_xvg_file("energy.xvg")
					E = pot[-10:].mean()/tot.max()
					if E < energy_limit:
						condition = True
					else:
						condition = False

					if condition: 
						doneflag[k] = True
						log = f"Run {config_name} from {self.name} did meet energy criteria. {E} > {energy_limit} and is finished ({counts[k]}) at {datetime.datetime.now()}.\n"
						with open(log_ext_path,'a') as f:
							f.write(log)

						
					if not condition:
						p = self.exp_sim2.gmx_grunpp(configuration,config_name,ionize=True,par=True,topol_file="../topol.top")
						proc.append(p)
						counts[k] +=1
						log = f"Run {config_name} from {self.name} did not meet energy criteria. {E} > {energy_limit} and was extended ({counts[k]}).\n"
						with open(log_path,'a') as f:
							f.write(log)
						with open(log_ext_path,'a') as f:
							f.write(log)

	
				for p in proc: p.wait()

			os.system("rm \#mdout*")


	def export_trr(self,all_results):

		## Make directory for trr files
		os.chdir(self.result_path)
		trr_dir = f"{self.result_path}/trr_dir"
		os.mkdir("trr_dir")

		# Move strucutre
		os.chdir(self.simulation_path)
		shutil.copyfile(f"{self.name}.gro",f"{trr_dir}/structure.gro")

		# Go through directories and copy trr files 
		dirs = [x for x in sorted(os.listdir(self.simulation_path)) if len(x.split('.')) == 1]
		for dir in dirs:
			os.chdir(self.simulation_path)
			os.chdir(dir)

			try: # Look for trr file
				files = os.listdir(os.getcwd())
				file = [x for x in files if x.split('.')[-1] == "trr"][0]
				shutil.copyfile(file,f"{trr_dir}/{file}")
			except: # If no trr file is found
				log = f"Can not find trr file for {self.name} for run in {dir}, simulation most likely crashed or never started.\n"
				with open(log_path,'a') as f:
					f.write(log)
				with open(log_ext_path,'a') as f:
					f.write(log)


		configurations = sorted([x for x in os.listdir(trr_dir) if x.split('.')[-1] == "trr"])

		# TODO: change output format.
		## load trr files to extract velocity vectors 
		data = []
		for l, configuration in enumerate(configurations):  
			U = md.Universe(f"{trr_dir}/structure.gro",f"{trr_dir}/{configuration}")
			ag = U.atoms.select_atoms("all")
			idx = [ag[k].index for k in range(len(ag))]
        
			velocity = U.trajectory[-1].velocities[idx]
			velocity = [(x/np.linalg.norm(x)) for x in velocity] 
			velocity = np.array(velocity, dtype=np.double)

			data.append(velocity)
		data = np.array(data)

		pkl_name = f"{self.name}_all"
		pkl_file = f"{self.result_path}/{pkl_name}"

		with open(pkl_file,'wb') as fp:
			pickle.dump(data,fp)

		os.chdir(trr_dir)
		os.system("rm .*")

	def cleanup_files(self):
		os.chdir(self.result_path)
		os.system("rm -r *")

	def collect_output(self):
		trr_dir = f"{self.result_path}/trr_dir"

		# Open H5 file
		with h5py.File(h5_path,"a") as file:
			group = file[f"{params.num_photons:.0E}/{model.name}/pdb_files"]

 			# Go through directories and extract last frame from trr files 
			dirs = [x for x in sorted(os.listdir(self.simulation_path)) if len(x.split('.')) == 1]
			for k,dir in enumerate(dirs):
				os.chdir(self.simulation_path)
				os.chdir(dir)

				# Find the trr and tpr file
				files = os.listdir(os.getcwd())
				trr_file = [x for x in files if x.split('.')[-1] == "trr"][0]
				tpr_file = [x for x in files if x.split('.')[-1] == "tpr"][0]

				# extract the last frame as a pdb
				t_end = self.exp_sim2.ts*self.exp_sim2.nsteps
				pdb_file = f"{trr_file.split('.')[0]}.pdb"
				cmd = f"echo '0' | {self.conf_sim.trjconv} -f {trr_file} -s {tpr_file} -b {t_end -5*self.exp_sim2.ts} -dump {t_end} -o ./{pdb_file}"  
				os.system(cmd)


				try: # If a simulation crashes for some reason this is where that usually will manifest.
					with open(pdb_file, "rb") as fp:
						pdb_data = fp.read()
					group.create_dataset(f"{k:05d}",data=pdb_data)

				except:
					with open(log_path,'a') as f:
						f.write(f"Last frame cannot be read for {self.name} on run {k:05d}, will use last existing frame.\n")


				
			configurations = sorted([x for x in os.listdir(trr_dir) if x.split('.')[-1] == "trr"])
			## load trr files to extract velocity vectors 
			data = []
			for l, configuration in enumerate(configurations):  
				U = md.Universe(f"{trr_dir}/structure.gro",f"{trr_dir}/{configuration}")
				ag = U.atoms.select_atoms("all")
				idx = [ag[k].index for k in range(len(ag))]
        
				velocity = U.trajectory[-1].velocities[idx]
				velocity = [(x/np.linalg.norm(x)) for x in velocity] 
				velocity = np.array(velocity, dtype=np.double)

				data.append(velocity)
			data = np.array(data)

			group = file[f"{params.num_photons:.0E}/{model.name}"]
			group.create_dataset("directions",data=data)











##############################################################
'''
HELPER FUNCTIONS
'''

def initialize_models(input_path,result_path,configuration_sim,explosion_sim1,explosion_sim2):
	inputs = os.listdir(input_path)

	models = []
	for item in inputs:
		name = item.split(".")[0]
		filepath = f"{input_path}/{item}"
		models.append(Model(name,filepath,result_path,configuration_sim,explosion_sim1,explosion_sim2))

	return models


def read_xvg_file(file_path):
    # Read the data from the file, skipping lines starting with '#' and '@'
    data = np.genfromtxt(file_path, comments='@',skip_header=10)

    # Separate the columns into NumPy arrays
    time = data[:,0]
    pot  = data[:,1]
    tot  = data[:,2]

    return time, pot, tot

def generate_atomic_parameters(photon_energy,pdb_name,src,dst):
	os.system(f"python3 {script_path}/generate_atomic_parameters.py {photon_energy} {pdb_name} {src} {dst}")


'''
______________________________________________________________
 ___      ___      __       __   _____  ___   
|"  \    /"  |    /""\     |" \ (\"   \|"  \  
 \   \  //   |   /    \    ||  ||.\\   \    | 
 /\\  \/.    |  /' /\  \   |:  ||: \.   \\  | 
|: \.        | //  __'  \  |.  ||.  \    \. | 
|.  \    /:  |/   /  \\  \ /\  |\    \    \ | 
|___|\__/|___(___/    \___|__\_|_)___|\____\) 
                                             


'''


# Path setup 
if not 1 < len(sys.argv) <= 3:
	raise Exception(f"Invalid number of arguments. {len(sys.argv)} given, expected 1 or 2.\n Usage: python3 exploder.py 'run_name' 'parameter_file'(optional)") 

run_name = sys.argv[1]
run_path = os.getcwd()
main_path  = os.path.realpath(os.path.dirname(__file__))

input_path 	= os.path.join(main_path,"input")
script_path = os.path.join(main_path,"scripts")
mdp_path 		= os.path.join(main_path,"mdp_files")
mdp_conf 		= os.path.join(mdp_path,"conf.mdp")
mdp_exp1  	= os.path.join(mdp_path,"exp1.mdp")
mdp_exp2  	= os.path.join(mdp_path,"exp2.mdp")




# TODO: better parameters handling
if len(sys.argv) == 3: # parameter file given
	parameter_file = str(sys.argv[2])
	log = f"Using provided parameter file: {parameter_file}.py"
	print(log)

else: # using default parameters
	log = f"No parameter file given, using 'params.py' in {main_path} as default"
	print(log)

	sys.path.append(main_path)
	parameter_file = "params"

params = __import__(parameter_file)

if params.result_dir_name not in os.listdir(params.result_path):
	os.chdir(params.result_path)
	os.mkdir(params.result_dir_name)


result_path = os.path.join(params.result_path,params.result_dir_name)
run_name_base = run_name
count = 1
while run_name in os.listdir(result_path):
	print(f"Runname '{run_name}' already exists, changing to {run_name_base}_{count}.")
	run_name = f"{run_name_base}_{count}"
	count = count + 1	

result_path = os.path.join(result_path,run_name)
os.mkdir(result_path)

# Save parameters
shutil.copyfile(f"{main_path}/{parameter_file}.py",f"{result_path}/parameters.txt")

### LOAD PARAMETERS FROM PARAMTER FILE ###

gmx_path 		 		 			= params.gmx_path
gmxp_path 		 		 		= params.gmxp_path # EXPLOSION PATH
extracted_frames 		  = params.extracted_frames
cleanup 							= params.cleanup
energy_limit					= params.energy_limit


# initialize Simulations
configuration_simulation = Simulation(
								path_to_mdp=mdp_conf,
								gmx_path=gmx_path,
								forcefield=params.forcefield,
								threads=params.threads_base,
								ts=params.ts_conf,
								nsteps=params.nsteps_conf,
								gaussian_peak=0, 			# FEL parameters are not used in config simulation
								num_photons=0,								
								sigma=0,
								focal_diameter=0,
								photon_energy=0,
								do_ff=0,
								do_ct=0,
								do_debye=0,
								do_log=0,
								read_states=0,
								double=False,
								MD_log=1,
								)

# This is the main simulation
explosion_simulation1 = Simulation(
								path_to_mdp=mdp_exp1,
								gmx_path=gmxp_path,
								forcefield=params.forcefield,
								threads=params.threads_exp,
								ts=params.ts_exp1,
								nsteps=params.nsteps_exp1,
								gaussian_peak=params.gaussian_peak,
								num_photons=params.num_photons,
								sigma=params.sigma,
								focal_diameter=params.focal_diameter,
								photon_energy=params.photon_energy,
								do_ff=params.do_ff1,
								do_ct=params.do_ct1,
								do_debye=params.do_debye1,								
								do_log=params.do_log1,
								read_states=params.read_states1,
								double=False,
								MD_log=params.MD_log,
								)

# Extension of simulation if needed
explosion_simulation2 = Simulation(
								path_to_mdp=mdp_exp2,
								gmx_path=gmxp_path,
								forcefield=params.forcefield,
								threads=params.threads_exp,
								ts=params.ts_exp2,
								nsteps=params.nsteps_exp2,
								gaussian_peak=0, # All FEL parameters 0, no new lazorz
								num_photons=0,
								sigma=0,
								focal_diameter=0,
								photon_energy=params.photon_energy,
								do_ff=params.do_ff2,
								do_ct=params.do_ct2,
								do_debye=params.do_debye2,			
								do_log=params.do_log2,
								read_states=params.read_states2,
								double=False,
								MD_log=params.MD_log,
								)

# Initialize Model(s)
models = initialize_models(input_path,
						result_path,
						configuration_simulation,
						explosion_simulation1,
						explosion_simulation2,
						)

# Initialize log file
log_path = f"{result_path}/{run_name}.log"
log_ext_path = f"{result_path}/{run_name}_extended.log"
log1 = f"Starting simulation on {datetime.datetime.now()}.\n"
log2 = f"Will run {params.extracted_frames} simulations using {params.threads_exp} threads for the following systems:\n"
with open(log_path,'w') as f:
	f.write(log1)
	f.write(log2)
	for model in models: f.write(f"{model.name}\n")
with open(log_ext_path,'w') as f:
	f.write(log1)
	f.write(log2)
	for model in models: f.write(f"{model.name}\n")


## Initialize HDF5 output file
h5_path = f"{result_path}/{run_name}.h5"
with h5py.File(h5_path,'w') as file:
	intensity_group = file.create_group(f"{params.num_photons:.0E}")
	intensity_group.create_dataset("num_photons",data=params.num_photons) 
	intensity_group.attrs["description"] = "Number of photons bombarded on the system."

	for model in models:
		molecule_group = intensity_group.create_group(f"{model.name}")
		pdb_group = molecule_group.create_group("pdb_files")



# "Actual" code starts here
for model in models:
	model.configuration_simulation()
	model.explosion_simulation()
	model.export_trr(result_path)
	model.collect_output()
	if cleanup: model.cleanup_files()

with open(log_path,'a') as f:
	f.write(f"Simulations finished at {datetime.datetime.now()}.")
