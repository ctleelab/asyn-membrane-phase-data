# gmx grommp to take .gro and .top from DOPC 3480ps to DPPC

import util
from util import base_path
from pathlib import Path
import subprocess
import MDAnalysis as mda

# setting paths
sys = 2
frame = int(458000)
initial_dim = "8x8x25"
final_dim = "large"
lipid = "DPPC"
ensemble = "NVT"
strain = "0.15strain"


script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent / "mdps")

analysis = util.analysis_path
analysis_curv_folder = f"curvature_selection"
system_analysis_folder = (
    analysis
    / analysis_curv_folder
    / ensemble
    / "buckled"
    / f"system{sys}-{initial_dim}-{lipid}-{strain}-NVT"
)
system_analysis_folder.mkdir(exist_ok=True)
min_folder = system_analysis_folder / "min"
# system_analysis_time_folder = f"system{sys}-{initial_dim}-{lipid}-{frame}ps-{ensemble}"

minimization_number = 7.0

gro_file = system_analysis_folder / f"{final_dim}-compression-{frame}ps-{strain}.gro"
top_file = system_analysis_folder / f"{final_dim}-system.top"
tpr_file = min_folder / f"minimization{minimization_number}.tpr"
mdp_file_min1 = mdp_path / f"step{minimization_number}_minimization.mdp"
output_file_min1 = min_folder / f"minimization{minimization_number}"


# tpr for first minimization step
minimize1_tpr = f"gmx grompp -f {mdp_file_min1} -c {gro_file} -r {gro_file} -p {top_file} -o {tpr_file}"
subprocess.run(minimize1_tpr, shell=True, check=True)

# Run mdrun for the first minimization step
minimize1_run = f"gmx mdrun -deffnm {output_file_min1} -v -nt 64 -rdd 1.5"
subprocess.run(minimize1_run, shell=True, check=True)

# Run the second minimization step if there is no error
minimization_number += 0.1
mdp_file_min2 = mdp_path / f"step{minimization_number}_minimization.mdp"
gro_file_min1 = min_folder / f"minimization7.0.gro"
output_file_min2 = min_folder / f"minimization{minimization_number}"
tpr_file_min2 = min_folder / f"minimization{minimization_number}.tpr"

# Create the .tpr file for the second minimization step
minimize2_tpr = f"gmx grompp -f {mdp_file_min2} -c {gro_file_min1} -r {gro_file_min1} -p {top_file} -o {tpr_file_min2}"
subprocess.run(minimize2_tpr, shell=True, check=True)

# Run mdrun for the second minimization step
minimize2_run = f"gmx mdrun -deffnm {output_file_min2} -v -nt 64"
subprocess.run(minimize2_run, shell=True, check=True)
