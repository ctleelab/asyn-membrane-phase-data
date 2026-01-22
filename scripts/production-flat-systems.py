#gmx grommp to take .gro and .top from DOPC 3480ps to DPPC 

import util
from util import base_path 
from pathlib import Path 
import subprocess
import MDAnalysis as mda

#inputs
sys = 6
initial_dim = "8x8x25"
final_dim = "medium"
lipid = "DOPC-DOPA"
ensemble = "NVT"
shape ="flat"

#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps")

analysis = util.analysis_path
analysis_curv_folder = f"curvature_selection"
system_analysis_folder = analysis/analysis_curv_folder/ensemble/shape/f"system{sys}-{initial_dim}-{lipid}-NVT"
equilibration_folder = system_analysis_folder / f"equil-{final_dim}"
equilibration_path = system_analysis_folder/ equilibration_folder

top_file = system_analysis_folder/ f"{final_dim}-system.top"
    
#create folder and path for production
production_folder = system_analysis_folder / f"prod-{final_dim}"
production_folder.mkdir(exist_ok = True)
production_number = 7.7

gro = equilibration_folder/ f"equilibration7.6.gro"

# Define files 
prod_mdp = mdp_path / f"step{production_number}_production.mdp"
tpr_file = production_folder / f"production{production_number}.tpr"
output_file_prod = production_folder/ f"production{production_number}"
cpt_file = equilibration_path / f"equilibration7.6.cpt"


#create .tpr
prod_tpr = f" gmx grompp -p {top_file} -f {prod_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
subprocess.run(prod_tpr, shell = True, check = True)

#run production
prod_mdrun = f" gmx mdrun -nt 40 -update gpu -gpu_id 0  -pin off -v -pinstride 1 -nstlist 100 -deffnm {output_file_prod}"
subprocess.run(prod_mdrun, shell = True, check = True)
