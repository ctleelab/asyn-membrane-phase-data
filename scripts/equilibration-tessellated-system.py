#goal run the equlibration for each system 1-10

#goal loop through each system folder
#take the minimization6.1 files 

#import libraries
import subprocess
from pathlib import Path 
from util import base_path


#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps") 


systems = [1,2,3,4,5,6]
size = "large"
for sys in systems:
    system_folder = f"system{sys}-8x8x25"
    system_path = base_path / system_folder
        
    #create folder and path for equilibration
    equilibration_folder = system_path / f"equil-{size}"
    equilibration_folder.mkdir(exist_ok = True)
    equilibration_path = system_path / equilibration_folder
    equilibration_number = 6.2
    previous_equilibration_number = round(equilibration_number - 0.1,1) 
    gro = system_path/ f"min-{size}"/f"minimization6.1.gro"
    #gro = equilibration_path/f"equilibration{previous_equilibration_number}.gro" #set if only running extentensions of equil 
    tesselation_folder = system_path/ "tessellation"

    while equilibration_number <= 6.9:
        # Define files 
        system_top = tesselation_folder/ f"{size}-system.top"
        equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
        tpr_file = equilibration_path / f"equilibration{equilibration_number}.tpr"
        output_file_equil = equilibration_path / f"equilibration{equilibration_number}"
        if equilibration_number == 6.2: #running 6.2 equilibration step seperately from the other because it does not require a .cpt file
            #create .tpr
            equil_tpr = f" gmx grompp -p {system_top} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -maxwarn 1"
            subprocess.run(equil_tpr, shell = True, check = True)

            #run equilibration
            equil_mdrun = f" gmx mdrun -nt 32 -update gpu -gpu_id 0  -pin off -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
            subprocess.run(equil_mdrun, shell = True, check = True)

            #input the equilibration gro file as new gro
            gro = equilibration_path / f"equilibration{equilibration_number}"

            equilibration_number += .1

        #runs 6.3 - 6.6 equilibration steps
        if equilibration_number <= 6.7: 
            #re-define files 
            equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
            tpr_file = equilibration_path / f"equilibration{equilibration_number}.tpr"
            output_file_equil = equilibration_path / f"equilibration{equilibration_number}"
            previous_equilibration_number = round(equilibration_number - 0.1,1) 
            cpt_file = equilibration_path / f"equilibration{previous_equilibration_number}.cpt"

            #create .tpr
            equil_tpr = f" gmx grompp -p {system_top} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
            subprocess.run(equil_tpr, shell = True, check = True)

            #run equilibration
            equil_mdrun = f" gmx mdrun -nt 32 -update gpu -gpu_id 0  -pin off -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
            subprocess.run(equil_mdrun, shell = True, check = True)

            #update gro file 
            gro = equilibration_path/ f"equilibration{equilibration_number}"

            equilibration_number = round(equilibration_number + 0.1,1) 

        #these run the longer equilibration steps 
        if equilibration_number >= 6.8: 
            
            #re-define files 
            equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
            tpr_file = equilibration_path / f"equilibration{equilibration_number}.tpr"
            output_file_equil = equilibration_path / f"equilibration{equilibration_number}"
            previous_equilibration_number = round(equilibration_number - 0.1,1) 
            cpt_file = equilibration_path / f"equilibration{previous_equilibration_number}.cpt"

            #create .tpr
            equil_tpr = f" gmx grompp -p {system_top} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
            subprocess.run(equil_tpr, shell = True, check = True)

            #run equilibration
            equil_mdrun = f" gmx mdrun -nt 32 -bonded gpu -gpu_id 0  -pin on -pinoffset 80  -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
            subprocess.run(equil_mdrun, shell = True, check = True)

            #update gro file 
            gro = equilibration_path/ f"equilibration{equilibration_number}"

            equilibration_number = round(equilibration_number + 0.1,1) 

            
    
        







