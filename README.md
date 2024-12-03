# **Description**

This repository allows for launching Molecular Dynamics (MD) simulations of HLA antigens from PDB structures. However, the current code could also be used to launch MD simulations of any type of structure that involves only protein structures.

# **Usage**

python pipeline_run_MD_HLAs.py <list_pdbs.txt>

Some warnings to bear in mind :

- The script utilizes the directory of the PDB file as working directory, so output will be generated in that directory. 

- This script assumes that the MD calculation will be launched in an OAR-based computation grid. If another resource/task manager
is used, please adapt the run_VMD function accordingly.

- This script needs the VMD-python package (version 3.1.4 tested).

- This script needs the NAMD3 executable and CHARMM36 topology and parameter files. Please change the variable "namd_dir" to 
the right directory pointing to these files (all these files must be in the same directory).

- Format example for 'list_pdbs.txt' (1 structure per line):
    /path/to/A_0201.pdb
    /path/to/A_1101.pdb
    /path/to/A_0202.pdb

# **Requirements**

- Python 3 (version 3.11.10 tested).
- VMD-python package (version 3.1.4 tested).
- NAMD3 (due to license issues, you must download a properly licensed NAMD3 from https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD and locate it at the *namd_dir* directory, version 3.0alpha13_multicore-CUDA tested).