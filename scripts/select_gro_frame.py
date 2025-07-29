#select a certain frame of trajectory from DOPC 50bar large 8x8x25
#then I want to change the DOPC for DPPC 
#then minimize, while holding the restraints of the phosphate of the lower leaflet'
#frame from vmd *100000 

import MDAnalysis as mda
import pathlib as Path 
import util
from util import base_path
import subprocess


systems = [1]
vmd_frame = 50
initial_dim = "8x8x25"
final_dim = 'large'
frame =  int(vmd_frame*1000*.02)

for sys in systems:
    system_folder = f"system{sys}-{initial_dim}"
    system_path = base_path/system_folder
    pressure_folder = "xzPcoupled-50bar-compression"
    pressure_path = system_path/pressure_folder
    file_name = f"{final_dim}-compression"
    input_name = f"{pressure_path}/{file_name}"
    frame_trr = f"echo '0'| gmx trjconv -f {input_name}.xtc -s {input_name}.tpr -o {input_name}-{frame}ps.gro -dump {frame}"
    subprocess.run(frame_trr, shell = True, check = True)


