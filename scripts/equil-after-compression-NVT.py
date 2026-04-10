# gmx grommp to take .gro and .top from DOPC 3480ps to DPPC

from util import base_path
from pathlib import Path
import subprocess

# inputs
final_dim = "8x8x25"
strains = ["strain.2", "flat"]
systems = [1,2,3,4,5,6]
size = "large"

# setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent / "mdps")
system_path = Path(script_dir.parent / "systems")
NVT_simulation_path = Path(script_dir.parent / "simulations-v2")


for system in systems:
    for strain in strains: 
        NVT_simulation_system_path = NVT_simulation_path/strain/f"system{system}"
        equilibration_folder = NVT_simulation_system_path/"equil"
        equilibration_folder.mkdir(exist_ok=True)
        top_file = system_path /f"system{system}-{final_dim}/tessellation/{size}-system.top"

        # set the starting equilibration number
        equilibration_number = 7.2
        gro = f"{ NVT_simulation_system_path}/min/minimization7.1.gro"

        while equilibration_number < 7.7:
            # Define files
            equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
            tpr_file = equilibration_folder / f"equilibration{equilibration_number}.tpr"
            output_file_equil = (
                equilibration_folder / f"equilibration{equilibration_number}"
            )

            if (
                equilibration_number == 7.2
            ):  # running 7.2 equilibration step seperately from the other because it does not require a .cpt file
                # create .tpr
                equil_tpr = f" gmx grompp -p {top_file} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -maxwarn 1"
                subprocess.run(equil_tpr, shell=True, check=True)

                # run equilibration
                equil_mdrun = f" gmx mdrun -nt 8 -update gpu -bonded gpu -gpu_id 0 -pin on -pinoffset 0 -v -pinstride 2 -nstlist 100 -deffnm {output_file_equil}"
                subprocess.run(equil_mdrun, shell=True, check=True)

                # input the equilibration gro file as new gro
                gro = equilibration_folder / f"equilibration{equilibration_number}.gro"

                equilibration_number += 0.1

            # runs 7.3 and above equilibration step. Seperate because you need a checkpoint file.
            if equilibration_number < 7.7:

                # re-define files
                equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
                tpr_file = equilibration_folder / f"equilibration{equilibration_number}.tpr"
                output_file_equil = (
                    equilibration_folder / f"equilibration{equilibration_number}"
                )
                previous_equilibration_number = round(equilibration_number - 0.1, 1)
                cpt_file = (
                    equilibration_folder
                    / f"equilibration{previous_equilibration_number}.cpt"
                )  # checkpoint from previous equilibrium run

                # creates .tpr to be used in the mdrun
                equil_tpr = f" gmx grompp -p {top_file} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
                subprocess.run(equil_tpr, shell=True, check=True)

                # run equilibration
                equil_mdrun = f" gmx mdrun -nt 8 -update gpu -bonded gpu -gpu_id 0 -pin on -pinoffset 0 -v -pinstride 2 -nstlist 100 -deffnm {output_file_equil}"
                subprocess.run(equil_mdrun, shell=True, check=True)

                # update gro file
                gro = equilibration_folder / f"equilibration{equilibration_number}.gro"

                equilibration_number = round(equilibration_number + 0.1, 1)
