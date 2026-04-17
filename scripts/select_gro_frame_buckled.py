import MDAnalysis as mda
import pathlib as Path 
import util
from util import base_path
import subprocess


initial_dim = "8x8x25"
final_dim = 'large'
configuration = "buckled"

#gives the frame number for each configurations
#buckled systems are at a strain of .2 
#flat systems are at the average Lx dimension 

buckled_systems = {1:"268740", 
                    2:"225420", 
                    3:"393180",
                    4:"367540",
                    5:"35100",
                    6:"190180" 
}

for system, frame in buckled_systems.items():
    configuration == "buckled"
    system_folder = f"system{system}-{initial_dim}"
    system_path = f"/scratch/local/casakurai/asyn-phase-binding-data/systems/{system_folder}"
    compression_folder = f"{system_path}/compression/xzPcoupled-35bar-compression-500ns-20psreadout"
    compression_file = f"{system_path}/compression/xzPcoupled-35bar-compression-500ns-20psreadout/large-compression"
    output_file = f"{compression_folder}/system{system}--buckled-strain.2-extracted"
    frame_trr = f"echo '0'| gmx trjconv -f {compression_file}.xtc -s {compression_file}.tpr -o {output_file}.gro -dump {frame}"
    subprocess.run(frame_trr, shell = True, check = True)