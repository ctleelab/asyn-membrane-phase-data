# gmx grommp to take .gro and .top from DOPC 3480ps to DPPC

import util
from util import base_path
from pathlib import Path
import subprocess
import MDAnalysis as mda

# setting paths
systems = [3,4,5,6]
initial_dim = "8x8x25"
final_dim = "large"
strain = "flat"


compression_3bar = [1,3,6]
compression_30bar = [5]
compression_35bar = [2]
compression_38bar = [4]


script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent / "mdps")
system_path = Path(script_dir.parent / "systems")
NVT_simulation_path = Path(script_dir.parent / "systems_NVT"/"simulations-v2")

for sys in systems:
    system_folder = f"{system_path}/system{sys}-{initial_dim}"

    compression_folder = f"{system_folder}/compression"
    tesselation_folder = f"{system_folder}/tessellation"

    if strain == "strain.2":
        if sys in compression_3bar:
            specific_compression_folder = f"{compression_folder}/xzPcoupled-3bar-compression-500ns-20psreadout"
        if sys in compression_30bar: 
            specific_compression_folder = f"{compression_folder}/xzPcoupled-30bar-compression-500ns-20psreadout"
        if sys in compression_35bar: 
            specific_compression_folder = f"{compression_folder}/xzPcoupled-35bar-compression-500ns-20psreadout"
        if sys in compression_38bar:
            specific_compression_folder = f"{compression_folder}/xzPcoupled-38bar-compression-500ns-20psreadout"
        gro_file = f"{specific_compression_folder}/system{sys}-buckled-strain.2-extracted.gro"

    if strain == "flat":
        gro_file = f"{system_folder}/equil-large/system{sys}-large-avg-flat-extractedgro.gro"

    NVT_analysis_folder = Path(
        f"{NVT_simulation_path}/{strain}/system{sys}"
    )
    NVT_analysis_folder.mkdir(exist_ok=True)
    min_folder = NVT_analysis_folder / "min"
    min_folder .mkdir(exist_ok=True)

    minimization_number = 7.0

    
    top_file = f"{tesselation_folder}/{final_dim}-system.top"
    tpr_file = f"{min_folder}/minimization{minimization_number}.tpr"
    mdp_file_min1 = f"{mdp_path}/step{minimization_number}_minimization.mdp"
    output_file_min1 = f"{min_folder}/minimization{minimization_number}"


    # tpr for first minimization step

    minimize1_tpr = f"gmx grompp -f {mdp_file_min1} -c {gro_file} -r {gro_file} -p {top_file} -o {tpr_file}"
    subprocess.run(minimize1_tpr, shell=True, check=True)

    # Run mdrun for the first minimization step
    minimize1_run = f"gmx mdrun -nt 8 -gpu_id 0 -pin on -pinoffset 0 -v  -pinstride 2 -rdd 1.5 -deffnm {output_file_min1} "
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
    minimize2_run = f"gmx mdrun -nt 8 -gpu_id 0 -pin on -pinoffset 0 -v -pinstride 2 -rdd 1.5 -deffnm {output_file_min2} "
    subprocess.run(minimize2_run, shell=True, check=True)
