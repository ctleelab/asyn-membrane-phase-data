# This script runs the production simulation for the buckled system. 
# It runs production step 7.8 and 7.9 which are both 600ns extensions of the 7.7 production run.

from pathlib import Path
import subprocess
import MDAnalysis as mda

# inputs
final_dim = "8x8x25"
strains = ["strain.2", "flat"]
systems = [1,2,3,4,5,6]
ensemble = "NVT"
size = "large"

# setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent / "mdps")
system_path = Path(script_dir.parent / "systems")
NVT_simulation_path = Path(script_dir.parent / "systems_NVT"/"simulations-v2")

for strain in strains: 
    for system in systems:
        NVT_simulation_system_path = NVT_simulation_path/strain/f"system{system}"
        equilibration_folder = NVT_simulation_system_path/"equil"
        equilibration_folder.mkdir(exist_ok=True)
        top_file = system_path /f"system{system}-{final_dim}/tessellation/{size}-system.top"


        production_folder = NVT_simulation_system_path / "prod"
    
        production_number = 7.8

        gro = production_folder  / f"production7.7.gro"

        # Define files
        prod_mdp = mdp_path / f"step{production_number}_production.mdp"
        tpr_file = production_folder / f"production{production_number}.tpr"
        output_file_equil = production_folder / f"production{production_number}"
        cpt_file = production_folder / f"production7.7.cpt"

        # create .tpr
        prod_tpr = f" gmx grompp -p {top_file} -f {prod_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
        subprocess.run(prod_tpr, shell=True, check=True)

        # run production
        prod_mdrun = f" gmx mdrun -nt 8 -update gpu -bonded gpu -gpu_id 0 -pin on -pinoffset 0 -v -pinstride 2 -nstlist 100 -deffnm {output_file_equil}"
        subprocess.run(prod_mdrun, shell=True, check=True)

        #next production run 
        production_number = 7.9
        
        gro = production_folder  / f"production7.8.gro"

        # Define files
        prod_mdp = mdp_path / f"step{production_number}_production.mdp"
        tpr_file = production_folder / f"production{production_number}.tpr"
        output_file_equil = production_folder / f"production{production_number}"
        cpt_file = production_folder / f"production7.8.cpt"

        # create .tpr
        prod_tpr = f" gmx grompp -p {top_file} -f {prod_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
        subprocess.run(prod_tpr, shell=True, check=True)

        # run production
        prod_mdrun = f" gmx mdrun -nt 8 -update gpu -bonded gpu -gpu_id 0 -pin on -pinoffset 0 -v -pinstride 2 -nstlist 100 -deffnm {output_file_equil}"
        subprocess.run(prod_mdrun, shell=True, check=True)