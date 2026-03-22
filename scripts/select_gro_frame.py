#select a certain frame of trajectory from DOPC 50bar large 8x8x25
#then I want to change the DOPC for DPPC 
#then minimize, while holding the restraints of the phosphate of the lower leaflet'
#frame from vmd *100000 

import MDAnalysis as mda
import pathlib as Path 
import util
from util import base_path
import subprocess


systems = [4]
initial_dim = "8x8x25"
final_dim = 'large'
configuration = "buckled"


#use if you're pulling from a particular time (ps)
frame =  int(270000)

if configuration == "flat":
    for sys in systems:
        system_folder = f"system{sys}-{initial_dim}"
        system_path = base_path/system_folder
        equilibration_folder = f"equil-{final_dim}"
        equilibration_path = system_path/ equilibration_folder
        output_file = f"system{systems}-{final_dim}-flat-extractedgro-{frame}ps"
        file_name = f"equilibration6.7"
        frame_trr = f"echo '0'| gmx trjconv -f {equilibration_path}/{file_name}.xtc -s {equilibration_path}/{file_name}.tpr -o {equilibration_path}/{output_file}.gro -dump {frame}"
        subprocess.run(frame_trr, shell = True, check = True)

if configuration == "buckled":
    for sys in systems:
        system_folder = f"system{sys}"
        system_path = f"scratch/local/casakurai/asyn-phase-binding-data/compression-data/{system_folder}"
        compression_file = system_path/f"{final_dim}-compression"
        output_file = f"system{systems}-{final_dim}-buckled-extractedgro-{frame}ps"
        file_name = f"equilibration6.7"
        frame_trr = f"echo '0'| gmx trjconv -f {compression_file}.xtc -s {compression_file}.tpr -o {compression_file}.gro -dump {frame}"
        subprocess.run(frame_trr, shell = True, check = True)



