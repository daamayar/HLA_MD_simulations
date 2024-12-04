#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script run MD simulations of HLA antigens.

USAGE: python pipeline_run_MD_HLAs.py <list_pdbs.txt>

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

Author: {0}
Email: {1}


'''


import os
import sys
from vmd import evaltcl
from vmd import molecule
import re
import multiprocessing as mp
import traceback
from vmd import atomsel


"""
PLEASE DEFINE THE "namd_dir" VARIABLE JUST BELLOW:
"""
namd_dir = os.getcwd() + '/namd_dir' # os.getcwd()

CWD = os.getcwd()

__author__ = "Diego Amaya"
__email__ = "diego.amaya-ramirez@inria.fr, diaamayaram@unal.edu.co"

USAGE = __doc__.format(__author__, __email__)


def check_input(args):
    try:
        """This function checks if to read from stdin/file and validates user input."""
        if not len(args):
            # Read from pipe
            if not sys.stdin.isatty():
                input_file = sys.stdin
            else:
                raise Exception
        elif len(args) == 1:
            if not os.path.isfile(args[0]):
                raise Exception
            input_file = open(args[0], 'r')
        else:
            raise Exception
    except:
        sys.exit(1)
    return input_file


def run_VMD(pdb_file, **kws):
    """
    pdbid = molecule.load('pdb', f'{pdb_file}')
    evaltcl('package require psfgen')
    evaltcl('resetpsf')
    evaltcl(f'topology {namd_dir}/top_all36_prot.rtf')

    #evaltcl(f"autopsf -mol {pdbid} -protein -prefix {pdb_file.split('.')[0]}")
    evaltcl("guesscoord")
    #evaltcl(f"writepdb {pdb_file.split('.')[0]}_autopsf.pdb")
    #evaltcl(f"writepsf {pdb_file.split('.')[0]}_autopsf.psf")
    molecule.write(pdbid, "pdb", f"{pdb_file.split('.')[0]}_autopsf.pdb")
    molecule.write(pdbid, "psf", f"{pdb_file.split('.')[0]}_autopsf.psf")
    molecule.delete(pdbid)
    evaltcl('resetpsf')
    pdbid = molecule.load('pdb', f"{pdb_file.split('.')[0]}_autopsf.pdb", 'psf', f"{pdb_file.split('.')[0]}_autopsf.psf")
    
    evaltcl('package require solvate')
    evaltcl('package require autoionize')
    #liste = [f'{pdb_file.split('.')[0]}']
    structure = f"{pdb_file.split('.')[0]}"

    #for structure in liste:
    os.system(f'mv {structure}_autopsf.psf {structure}.psf')
    os.system(f'mv {structure}_autopsf.pdb {structure}.pdb')
    molecule.delete(pdbid)

    evaltcl('resetpsf')
    pdbid = molecule.load('pdb', f'{structure}.pdb', 'psf', f'{structure}.psf')
    evaltcl(f'set all [atomselect {pdbid} all]')
    evaltcl(f'$all moveby [vecinvert [measure center $all weight mass]]')
    evaltcl('resetpsf')
    evaltcl(f'solvate {structure}.psf {structure}.pdb -t 20 -o {structure}_wb')
    evaltcl(f'autoionize -psf {structure}_wb.psf -pdb {structure}_wb.pdb -neutralize -o {structure}_ionized ')
    molecule.delete(pdbid)
    """
    #evaltcl('resetpsf')
    structure = f"{pdb_file.split('.')[0]}"
    evaltcl(f'set pdb_file "{structure}"')
    evaltcl(f'set namd_dir "{namd_dir}"')
    evaltcl(f"source {CWD}/prep_MD.tcl")
    pdbid = molecule.load('pdb', f'{structure}_ionized.pdb', 'psf', f'{structure}_ionized.psf')
    evaltcl(f'set all [atomselect {pdbid} all]')
    cell_size = evaltcl('measure minmax $all')
    cell_size = re.findall('\{[^}]*\}', cell_size)
    cell_size = [elem.strip('{}') for elem in cell_size]
    cell_size = [elem.split(' ') for elem in cell_size]
    X_size = abs(float(cell_size[0][0])) + abs(float(cell_size[1][0]))
    Y_size = abs(float(cell_size[0][1])) + abs(float(cell_size[1][1]))
    Z_size = abs(float(cell_size[0][2])) + abs(float(cell_size[1][2]))
    
    sel = atomsel("all")
    sel.beta = 0.0
    sel = atomsel("protein")
    sel.beta = 1.0
    molecule.write(pdbid, "pdb", f"{structure}_ionized_allConstraints.pdb")
    molecule.delete(pdbid)

    evaltcl('resetpsf')
    pdbid = molecule.load('pdb', f'{structure}_ionized.pdb', 'psf', f'{structure}_ionized.psf')
    sel = atomsel("all")
    sel.beta = 0.0
    sel = atomsel("backbone")
    sel.beta = 1.0
    molecule.write(pdbid, "pdb", f"{structure}_ionized_bbConstraints.pdb")
    molecule.delete(pdbid)

    evaltcl('resetpsf')
    pdbid = molecule.load('pdb', f'{structure}_ionized.pdb', 'psf', f'{structure}_ionized.psf')
    min_conf = open(f'{namd_dir}/general_min.conf', 'r')
    new_min_conf = open(f'{structure}_min.conf', 'w')
    for line in min_conf:
        if '<file_name>' in line:
            new_min_conf.write(line.replace('<file_name>', f'{structure}_ionized'))
        elif '<X_size>' in line:
            new_min_conf.write(line.replace('<X_size>', f'{X_size}'))
        elif '<Y_size>' in line:
            new_min_conf.write(line.replace('<Y_size>', f'{Y_size}'))
        elif '<Z_size>' in line:
            new_min_conf.write(line.replace('<Z_size>', f'{Z_size}'))
        elif '<center_X>   <center_Y>   <center_Z>' in line:
            new_min_conf.write(line.replace('<center_X>   <center_Y>   <center_Z>',
                    f'{(float(cell_size[0][0]) + float(cell_size[1][0]))/2}   {(float(cell_size[0][1]) + float(cell_size[1][1]))/2}   {(float(cell_size[0][2]) + float(cell_size[1][2]))/2}'))
        elif '<namd_dir>' in line:
            new_min_conf.write(line.replace('<namd_dir>', f'{namd_dir}'))
        else:
            new_min_conf.write(line)
    
    min_conf.close()
    new_min_conf.close()

    allConstrains_conf = open(f'{namd_dir}/general_allConstrains.conf', 'r')
    new_allConstrains_conf = open(f'{structure}_allConstrains.conf', 'w')
    for line in allConstrains_conf:
        if '<file_name>' in line:
            new_allConstrains_conf.write(line.replace('<file_name>', f'{structure}_ionized'))
        elif '<namd_dir>' in line:
            new_allConstrains_conf.write(line.replace('<namd_dir>', f'{namd_dir}'))
        else:
            new_allConstrains_conf.write(line)
    allConstrains_conf.close()
    new_allConstrains_conf.close()

    bbConstrains_conf = open(f'{namd_dir}/general_bbConstrains.conf', 'r')
    new_bbConstrains_conf = open(f'{structure}_bbConstrains.conf', 'w')
    for line in bbConstrains_conf:
        if '<file_name>' in line:
            new_bbConstrains_conf.write(line.replace('<file_name>', f'{structure}_ionized'))
        elif '<namd_dir>' in line:
            new_bbConstrains_conf.write(line.replace('<namd_dir>', f'{namd_dir}'))
        else:
            new_bbConstrains_conf.write(line)
    bbConstrains_conf.close()
    new_bbConstrains_conf.close()

    production_conf = open(f'{namd_dir}/general_production.conf', 'r')
    new_production_conf = open(f'{structure}_production.conf', 'w')
    for line in production_conf:
        if '<file_name>' in line:
            new_production_conf.write(line.replace('<file_name>', f'{structure}_ionized'))
        elif '<namd_dir>' in line:
            new_production_conf.write(line.replace('<namd_dir>', f'{namd_dir}'))
        else:
            new_production_conf.write(line)
    production_conf.close()
    new_production_conf.close()

    os.system('chmod +x *.conf')
    run = open(f'run_MD_{structure}.sh', 'w')
    run.writelines([
            '#! /bin/bash\n',
            f'bash -l -c "{namd_dir}/namd3 +idlepoll +p {int(mp.cpu_count()/4)} +setcpuaffinity +devices 0 {structure}_min.conf > {structure}_min.log"\n',
            f'bash -l -c "{namd_dir}/namd3 +idlepoll +p 1 +setcpuaffinity +devices 0 {structure}_allConstrains.conf > {structure}_allConstrains.log"\n',
            f'bash -l -c "{namd_dir}/namd3 +idlepoll +p 1 +setcpuaffinity +devices 0 {structure}_bbConstrains.conf > {structure}_bbConstrains.log"\n',
            f'bash -l -c "{namd_dir}/namd3 +idlepoll +p 1 +setcpuaffinity +devices 0 {structure}_production.conf > {structure}_production.log"'
            ])
    run.close()
    os.system(f'chmod +x run_MD_{structure}.sh')
    cluster_list = ['graffiti' , 'gruss', 'grue', 'grat']
    os.system(f"oarsub --notify mail:diego.amaya-ramirez@inria.fr -n {pdb_file.split('.')[0]} -q production -p 'cluster in ({', '.join(cluster_list)})' -l host=1/gpu=1,walltime=96 ./run_MD_{structure}.sh")
    molecule.delete(pdbid)


if __name__ == "__main__":
    try:
        input_file = check_input(sys.argv[1:])
        
        for line in input_file:
            path_input, file = os.path.split(os.path.abspath(line))
            os.chdir(f'{path_input}')
            run_VMD(file.strip())
        
        input_file.close()
    except:
        sys.stderr.write(USAGE)
        traceback.print_exc()
        sys.exit(1)
