# gmx grommp to take .gro and .top from DOPC 3480ps to DPPC

import util
from util import base_path
from pathlib import Path
import subprocess
import MDAnalysis as mda

# inputs
final_dim = "8x8x25"
strain = "0.15strain"
systems = [
    f"system1-{final_dim}-DOPC-{strain}-NVT",
    f"system2-{final_dim}-DPPC-{strain}-NVT",
    f"system3-{final_dim}-DOPC-DOPS-{strain}-NVT",
    f"system4-{final_dim}-DPPC-DPPS-{strain}-NVT",
    f"system5-{final_dim}-DPPC-DPPA-{strain}-NVT",
    f"system6-{final_dim}-DOPC-DOPA-{strain}-NVT",
]
ensemble = "NVT"
size = "large"

# setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent / "mdps")

analysis = util.analysis_path
analysis_curv_folder = f"curvature_selection"

for system in systems:
    system_analysis_folder = (
        analysis / analysis_curv_folder / ensemble / "buckled" / system
    )
    equilibration_folder = system_analysis_folder / "equil"
    top_file = system_analysis_folder / f"{size}-system.top"

    # create folder and path for production
    production_folder = system_analysis_folder / "prod"
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
    equil_mdrun = f" gmx mdrun -nt 64  -update gpu -gpu_id 0  -pin off -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
    subprocess.run(equil_mdrun, shell=True, check=True)
