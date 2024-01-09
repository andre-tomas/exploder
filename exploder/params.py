# DEFAULT PARAMETER FILE
# This file must be in same location as script is run from!

# gmx paths
gmx_path = '/home/tomas/programs/gromacs-4.5.4-unmod/bin' 	# Gromacs Unmodified (4.5.4) bin directory
gmxp_path =	'/home/tomas/programs/gromacs-mc/bin'# Gromacs Explosion  (4.5.4) bin directory

# Results paths
result_dir_name = "results"
result_path ="/home/tomas/Documents/jobb/bra_att_ha_koder/exploder" # Where result directory will be saved
cleanup = False # if True keep only reduced results (remove all big trr files)
# The collected 'reduced' results are put in .../result_path/output


# Forcefield parameters
forcefield ="charmm36-mar2019-Fe-S"  # Forcefield for HiPIP
#forcefield = "CHARMM"
water = "spce"		  		# Water model used (if needed)

# Number of simulations
threads_base = 1
threads_exp  = 1	# Sometimes too many threads crashes simulation for exp sims
extracted_frames = 1 # Number of starting configurations extracted per system
MD_log = 100 # How often to write to the trr files (in steps)


# Configuration sims 
nsteps_conf = 1_000 # 1_000_000 # should be atleast 10*extracted frames
ts_conf     = 0.001 

# Explosion sim 1
nsteps_exp1 = 50_000
ts_exp1     = 0.000_001  # 0.001 = 1.0fs, 1 = 1ps, 0.05fs = 50as = 0.000_05 ps, 0.000_001 ps = 1as
do_ff1		=1   # alter forcefield userint1
do_cg1		=1   # Do charge transfer  userint2
do_debye1	=1   # Use debye shielding  userint3
do_log1		=1  # Write log files (Use for debugging, big performance drop) userint5
read_states1=0   # (0) Start from ground states, (1) Read states and charges from file. userint4

# FEL parameters
guassian_peak = 0.0 # [ps] userreal1 
num_photons   = 1e11 # userreal2
FWHM          = 0.010616 # userreal3
focal_diameter= 100 # [nm] userreal4


# Explosion sim 2
nsteps_exp2 = 10_000
ts_exp2     = 0.000_001
do_ff2		=1 	  # alter forcefield
do_cg2		=1    # Do charge transfer 
do_debye2	=1    # Use debye shielding 
do_log2		=1    # Write log files (Use for debugging, big performance drop)
read_states2=1    # (0) Start from ground states, (1) Read states and charges from file.
energy_limit = 0.90
# No FEL parameters for second explosion.
# This simulation runs as a contiunuation if the previous one did not fulfill the energy condition.
# E_pot/E_tot < energy_limit


#####
# How to start:
## 1. conda activiate HiPIP
## 2. nohup nice -n 19 python3 exploder.py "runname" params > output.log 2>&1 &
## 3. (optional)disown -h "process idq"
