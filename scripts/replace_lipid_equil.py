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

#use if you're pulling from a particular time (ps)
frame =  int(85620)
initial_dim = "8x8x25"
final_dim = "large"
compression = "3bar"
lipid = "DPPC"
sim_time = 200
ensemble = "NVT"

#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps")

analysis = util.analysis_path
analysis_curv_folder = f"curvature_selection"
system_analysis_folder = analysis/analysis_curv_folder/f"{compression}_{sim_time}ns"/f"system{sys}-{initial_dim}-{lipid}-{frame}ps-NVT"
system_analysis_folder.mkdir(exist_ok = True)
system_analysis_time_folder = f"system{sys}-{initial_dim}-{lipid}-{frame}ps-{ensemble}"
analysis_curv_dir = analysis/analysis_curv_folder/f"{compression}_{sim_time}ns"/f"{system_analysis_folder}"/f"{system_analysis_time_folder}"
analysis_curv_dir.mkdir(exist_ok = True)


top_file = analysis_curv_dir / f"{final_dim}-system.top"
    
#create folder and path for equilibration
equilibration_folder = analysis_curv_dir / f"equil"
equilibration_folder.mkdir(exist_ok = True)
equilibration_path = analysis_curv_dir / equilibration_folder
equilibration_number = 7.2
gro = analysis_curv_dir / f"minimization7.1.gro"

while equilibration_number <= 7.7:
    # Define files 
    equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
    tpr_file = equilibration_path / f"equilibration{equilibration_number}.tpr"
    output_file_equil = equilibration_path / f"equilibration{equilibration_number}"
    if equilibration_number == 7.2: #running 7.2 equilibration step seperately from the other because it does not require a .cpt file
        #create .tpr
        equil_tpr = f" gmx grompp -p {top_file} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -maxwarn 1"
        subprocess.run(equil_tpr, shell = True, check = True)

        #run equilibration
        equil_mdrun = f" gmx mdrun -nt 4 -update gpu -gpu_id 0  -pin on -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
        subprocess.run(equil_mdrun, shell = True, check = True)

        #input the equilibration gro file as new gro
        gro = equilibration_path / f"equilibration{equilibration_number}"

        equilibration_number += .1
    #runs 6.3 and above equilibration step 
    if equilibration_number < 7.7: 
        #re-define files 
        equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
        tpr_file = equilibration_path / f"equilibration{equilibration_number}.tpr"
        output_file_equil = equilibration_path / f"equilibration{equilibration_number}"
        previous_equilibration_number = round(equilibration_number - 0.1,1) 
        cpt_file = equilibration_path / f"equilibration{previous_equilibration_number}.cpt"

        #create .tpr
        equil_tpr = f" gmx grompp -p {top_file} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
        subprocess.run(equil_tpr, shell = True, check = True)

        #run equilibration
        equil_mdrun = f" gmx mdrun -nt 4 -update gpu -gpu_id 0  -pin on -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
        subprocess.run(equil_mdrun, shell = True, check = True)

        #update gro file 
        gro = equilibration_path/ f"equilibration{equilibration_number}"

        equilibration_number = round(equilibration_number + 0.1,1) 

    # if equilibration_number == 7.7: 
    #     #re-define files 
    #     equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
    #     tpr_file = equilibration_path / f"equilibration{equilibration_number}.tpr"
    #     output_file_equil = equilibration_path / f"equilibration{equilibration_number}"
    #     previous_equilibration_number = round(equilibration_number - 0.1,1) 
    #     cpt_file = equilibration_path / f"equilibration{previous_equilibration_number}.cpt"

    #     #create .tpr
    #     equil_tpr = f" gmx grompp -p {top_file} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
    #     subprocess.run(equil_tpr, shell = True, check = True)

    #     #run equilibration
    #     equil_mdrun = f" gmx mdrun -nt 8 -update gpu -gpu_id 0  -pin on -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
    #     subprocess.run(equil_mdrun, shell = True, check = True)

    #     #update gro file 
    #     gro = equilibration_path/ f"equilibration{equilibration_number}"

    #     equilibration_number = round(equilibration_number + 0.1,1) 
        

    







