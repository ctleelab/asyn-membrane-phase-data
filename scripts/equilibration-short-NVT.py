#goal run the equlibration for each system 1-10

#goal loop through each system folder
#take the minimization6.1 files 

#import libraries
import subprocess
from pathlib import Path 
from util import base_path
from util import analysis_path


#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps") 


systems = 5
final_dim = "medium"
time = "220"
lipids = "DPPC-DPPA"

system_folder = f"system{systems}-8x8x25"
system_path = base_path / system_folder
    
#create folder and path for equilibration
equilibration_path = system_path / f"equil-{final_dim}"
tessellation_path = system_path / "tessellation"
gro_file = f"system{systems}-{final_dim}-flat-extractedgro-{time}ps.gro"
tesselation_folder = system_path/ "tessellation"
system_analysis_folder = analysis_path/"curvature_selection"/"NVT"/f"system{systems}-8x8x25-{lipids}-NVT"/"equil"

#creates a new .tpr file with DPPC instead of DOPC 
equilibration_number = 7.6

top_file = tessellation_path / f"{final_dim}-system.top"
tpr_file = equilibration_path / f"equilibration{equilibration_number}.tpr"
mdp_file_min1 = mdp_path / f"step{equilibration_number}_equilibration.mdp"
output_file_min1 = equilibration_path/ f"equilibration{equilibration_number}"


#tpr for first equil step 
equil1_tpr = f"gmx grompp -f {mdp_file_min1} -c {equilibration_path}/{gro_file} -r {equilibration_path}/{gro_file} -p {top_file} -o {tpr_file}"
subprocess.run(equil1_tpr, shell=True, check=True)

# Run mdrun for the equil
equil1_run = f"gmx mdrun -deffnm {output_file_min1} -v -nt 64 -rdd 1.5"
subprocess.run(equil1_run, shell=True, check=True)



        







