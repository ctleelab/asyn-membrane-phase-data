#load in the systems that has base paths 
import util as util
from util import base_path

import subprocess 

from pathlib import Path


#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps") 

#systems 
systems = [1,2]


for sys in systems: 
    # Reset the minimization number, needs to start at 6.0 for each system
    minimization_number = 6.0 
    # Define system folder and path    
    system_folder = f"system{sys}"
    system_path = base_path / system_folder
    
    # Define files 
    gro_file = system_path / "neutral.gro"
    top_file = system_path / "system.top"
    tpr_file = system_path / f"minimization{minimization_number}.tpr"
    mdp_file_min1 = mdp_path / f"step{minimization_number}_minimization.mdp"
    output_file_min1 = system_path / f"minimization{minimization_number}"

    try:
        # Create the .tpr file for the first minimization step
        minimize1_tpr = f"gmx grompp -f {mdp_file_min1} -c {gro_file} -r {gro_file} -p {top_file} -o {tpr_file}"
        subprocess.run(minimize1_tpr, shell=True, check=True)

        # Run mdrun for the first minimization step
        minimize1_run = f"gmx mdrun -deffnm {output_file_min1} -v -nt 64 -rdd 1.5"
        subprocess.run(minimize1_run, shell=True, check=True)

        # Run the second minimization step if there is no error
        minimization_number += 0.1
        mdp_file_min2 = mdp_path / f"step{minimization_number}_minimization.mdp"
        gro_file_min1 = system_path / f"minimization6.0.gro"
        output_file_min2 = system_path / f"minimization{minimization_number}"
        tpr_file_min2 = system_path / f"minimization{minimization_number}.tpr"

        # Create the .tpr file for the second minimization step
        minimize2_tpr = f"gmx grompp -f {mdp_file_min2} -c {gro_file_min1} -r {gro_file_min1} -p {top_file} -o {tpr_file_min2}"
        subprocess.run(minimize2_tpr, shell=True, check=True)

        # Run mdrun for the second minimization step
        minimize2_run = f"gmx mdrun -deffnm {output_file_min2} -v -nt 64"
        subprocess.run(minimize2_run, shell=True, check=True)

    except subprocess.CalledProcessError as e:
        print(f"Error occurred in system {sys} during minimization step {minimization_number}. Continuing to next system.")
        # Save the error to a text file
        error_file = system_path / "minimization_error.txt"
        with open(error_file, "w") as f:
            f.write(f"Error occurred in system {sys} during minimization step {minimization_number}.\n")


