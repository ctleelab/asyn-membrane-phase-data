#gmx grommp to take .gro and .top from DOPC 3480ps to DPPC 

import util
from util import base_path 
from pathlib import Path 
import subprocess
import MDAnalysis as mda

#inputs
sys = 1

#use if you are looking directly at the vmd simulation 
# vmd_frame = 50
# frame =  int(vmd_frame*1000*.02)

##use if you're pulling from a particular time (ps)
 
frame =  int(160140)
initial_dim = "8x8x25"
final_dim = "large"
lipid = "DOPC"
ensemble = "NVT"

#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps")

analysis = util.analysis_path
analysis_curv_folder = f"curvature_selection"
system_analysis_folder = analysis/analysis_curv_folder/ensemble/f"system{sys}-{initial_dim}-{lipid}-{frame}ps-NVT"
equil_system_folder = system_analysis_folder / "equil"
equilibration_folder = system_analysis_folder / f"equil"
equilibration_path = system_analysis_folder/ equilibration_folder

top_file = system_analysis_folder/ f"{final_dim}-system.top"
    
#create folder and path for production
production_folder = system_analysis_folder / "prod"
production_folder.mkdir(exist_ok = True)
production_number = 7.7

gro = equil_system_folder/ f"equilibration7.6.gro"

# Define files 
equil_mdp = mdp_path / f"step{production_number}_production.mdp"
tpr_file = production_folder / f"production{production_number}.tpr"
output_file_equil = production_folder/ f"production{production_number}"
cpt_file = equilibration_path / f"equilibration7.6.cpt"


#create .tpr
equil_tpr = f" gmx grompp -p {top_file} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
subprocess.run(equil_tpr, shell = True, check = True)

#run production
equil_mdrun = f" gmx mdrun -nt 4 -update gpu -gpu_id 0  -pin on -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
subprocess.run(equil_mdrun, shell = True, check = True)
