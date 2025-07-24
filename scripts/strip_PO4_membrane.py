#steps
#will want to average through out the length of the x-dimension, then how tight the variation in curvature
#first extract the phosphates from the trajectories 
#define upper and lower leaflet 
#create a function that represents one of the monolayers 
#then calculate local curvature


import MDAnalysis as mda
import matplotlib as plt
import numpy as np
import util
from util import base_path 

#library for showing progress bar for runs
from tqdm.auto import tqdm
#library that does parallel processing 
from tqdm.contrib.concurrent import process_map 
from pathlib import Path 
import os 
import subprocess 

script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps") 


initial_dim = "6x6x12"
final_dims = ['large']

#figures be created in current working directory, which should be scipts
figures = Path("figures")
figures.mkdir(parents=True, exist_ok=True)
figures_curvature = figures/"local_curvature"
figures_curvature.mkdir(parents=True, exist_ok=True)


def _strip_trajectory(job) -> None:
    sys, dim = job
    analysis = util.analysis_path
    analysis_curv_folder = f"local_curvature"
    analysis_curv_dir = analysis/analysis_curv_folder 
    analysis_curv_dir.mkdir(parents = True, exist_ok = True)

    os.chdir(analysis_curv_dir)
    system_folder = f"system{sys}-{initial_dim}"
    system_path = base_path/system_folder

    if sys == 1:
        pressure_folder = "100bar-xzPcoupled-compression"
    else:
        pressure_folder = "200bar-xzPcoupled-compression"
    pressure_path = system_path/pressure_folder

    for dim in final_dims:
        file_name = f"{dim}-compression"

        #create trajectory file with only PO4
        trjconv_cmd = f"echo '2' | gmx trjconv -f {pressure_path}/{file_name}.gro -o {sys}-{dim}-po4_only.gro -n {sys}-{dim}-po4_membrane.ndx -s {sys}-{dim}-analysis.tpr"
        subprocess.run(trjconv_cmd, shell=True, check=True)


        
        #center the trajectory which is important b/c membrane can move through the simulation box during the simulation
        trjconv_cmd = f"echo '1 2' | gmx trjconv -f {pressure_path}/{file_name}.xtc -center -o {sys}-{dim}-po4_all.xtc -n {sys}-{dim}-po4_membrane.ndx -s {sys}-{dim}-analysis.tpr"
        subprocess.run(trjconv_cmd, shell=True, check=True)





systems = [1, 2]
labels = {1: "DOPC", 2: "DPPC"}



jobs = []

for sys in systems:
    #file setup for where files will be pulled from 
    system_folder = f"system{sys}-{initial_dim}"
    system_path = base_path/system_folder
    analysis = util.analysis_path
    analysis.mkdir(parents=True, exist_ok=True)
    analysis_curv_folder = f"local_curvature"
    analysis_curv_dir = analysis/analysis_curv_folder 
    analysis_curv_dir.mkdir(parents = True, exist_ok = True)

    #takes into account the different pressure required to bend membrane if DOPC vs DPPC 
    if sys == 1:
        pressure_folder = "100bar-xzPcoupled-compression"
    else:
        pressure_folder = "200bar-xzPcoupled-compression"
    pressure_path = system_path/pressure_folder
    for dim in final_dims:
        file_name = f"{dim}-compression"
        input_name = f"{pressure_path}/{file_name}"
        u =  mda.Universe(f"{system_path}/tessilation/{dim}-system.top", 
        f"{pressure_path}/{file_name}.gro", 
        topology_format = "ITP")

        membrane = u.atoms.select_atoms(util.membrane_sel)
        PO4= u.atoms.select_atoms(util.po4_sel)
        all_atoms = u.atoms 

        #only extract a certain interval of frames

        #creating an ndx file 
        with mda.selections.gromacs.SelectionWriter(f"{analysis_curv_dir}/{sys}-{dim}-po4_membrane.ndx", mode = "w")as ndx:
            ndx.write(all_atoms,name = "system")
            ndx.write(PO4, name = "po4")
            ndx.write(membrane, name = "membrane")

        #extract .gro from 10,000ps 
        interval_tpr = f"echo '0' | gmx trjconv -f {pressure_path}/{file_name}.xtc -s {pressure_path}/{file_name}.tpr -o {analysis_curv_dir}/{sys}-{dim}-trimmed.gro -dump {10000}"
        subprocess.run(interval_tpr, shell=True, check = True)


        #generate the tpr files with the above categorization of atoms
        tpr1 = f" gmx grompp -p {system_path}/tessilation/{dim}-system.top -f {mdp_path}/step6.8_compression.mdp -n {analysis_curv_dir}/{sys}-{dim}-po4_membrane.ndx -c {analysis_curv_dir}/{sys}-{dim}-trimmed.gro -o {analysis_curv_dir}/{sys}-{dim}-analysis.tpr -maxwarn 1"
        subprocess.run(tpr1,shell=True,check=True)

        jobs.append((sys,dim))

#running 
r = process_map(_strip_trajectory, jobs, max_workers=12)

