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

    # set the starting equilibration number
    equilibration_number = 7.2
    gro = system_analysis_folder / "min" / "minimization7.1.gro"

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
            equil_mdrun = f" gmx mdrun -nt 32 -update gpu -gpu_id 0  -pin off -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
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
            equil_mdrun = f" gmx mdrun -nt 32 -update gpu -gpu_id 0  -pin off -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
            subprocess.run(equil_mdrun, shell=True, check=True)

            # update gro file
            gro = equilibration_folder / f"equilibration{equilibration_number}.gro"

            equilibration_number = round(equilibration_number + 0.1, 1)

    # if equilibration_number == 7.7:
    #     #re-define files
    #     equil_mdp = mdp_path / f"step{equilibration_number}_equilibration.mdp"
    #     tpr_file = equilibration_folder / f"equilibration{equilibration_number}.tpr"
    #     output_file_equil = equilibration_folder / f"equilibration{equilibration_number}"
    #     previous_equilibration_number = round(equilibration_number - 0.1,1)
    #     cpt_file = equilibration_folder / f"equilibration{previous_equilibration_number}.cpt"

    #     #create .tpr
    #     equil_tpr = f" gmx grompp -p {top_file} -f {equil_mdp} -c {gro} -r {gro} -o {tpr_file} -t {cpt_file} -maxwarn 1"
    #     subprocess.run(equil_tpr, shell = True, check = True)

    #     #run equilibration
    #     equil_mdrun = f" gmx mdrun -nt 8 -update gpu -gpu_id 0  -pin on -v -pinstride 1 -nstlist 100 -deffnm {output_file_equil}"
    #     subprocess.run(equil_mdrun, shell = True, check = True)

    #     #update gro file
    #     gro = equilibration_folder/ f"equilibration{equilibration_number}"

    #     equilibration_number = round(equilibration_number + 0.1,1)
