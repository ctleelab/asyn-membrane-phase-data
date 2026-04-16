# gmx grommp to take .gro and .top from DOPC 3480ps to DPPC

import util
from util import base_path
from pathlib import Path
import subprocess

# inputs
final_dim = "8x8x25"
strains = ["strain.2", "flat"]
systems = [1,2,3]
ensemble = "NVT"
size = "large"

# setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent / "mdps")
system_path = Path(script_dir.parent / "systems")
NVT_simulation_path = Path(script_dir.parent / "simulations-v2")

for strain in strains:
    for system in systems:
        NVT_simulation_system_path = NVT_simulation_path/strain/f"system{system}"
        equilibration_folder = NVT_simulation_system_path/"equil"
        equilibration_folder.mkdir(exist_ok=True)
        top_file = system_path /f"system{system}-{final_dim}/tessellation/{size}-system.top"



        # create folder and path for production
        production_folder = NVT_simulation_system_path / "prod"
        production_folder.mkdir(exist_ok=True)
        production_number = 7.7

        gro = equilibration_folder / f"equilibration7.6.gro"

        # Define files
        equil_mdp = mdp_path / f"step{production_number}_production.mdp"
        tpr_file = production_folder / f"production{production_number}.tpr"
        output_file_equil = production_folder / f"production{production_number}"
        cpt_file = equilibration_folder / f"equilibration7.6.cpt"

        # create .tpr
        equil_tpr = f" gmx grompp -p {top_file} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
        subprocess.run(equil_tpr, shell=True, check=True)

        # run production
        equil_mdrun = f" gmx mdrun -nt 8 -update gpu -bonded gpu -gpu_id 0 -pin on -pinoffset 0 -v -pinstride 2 -nstlist 100 -deffnm {output_file_equil}"
        subprocess.run(equil_mdrun, shell=True, check=True)


