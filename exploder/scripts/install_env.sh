# Command that creates and environment with all dependancies for code to run
# Requires a (Ana-/Mini) conda installation
conda create --name explode python=3.8 numpy matplotlib MDAnalysis subprocess32 h5py openbabel -c conda-forge
