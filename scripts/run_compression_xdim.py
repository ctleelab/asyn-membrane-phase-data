from pathlib import Path
import util as util
from util import base_path 
import subprocess
from box_size import avg_box_size

#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps") 
size = "large"
systems = [1,2,3,4,5,6]
#dimensions = {'xsmall': "1 1 1", 'small': "2 1 1", 'medium': "3 1 1", 'large': "4 1 1"}
compression_3bar = [1,3,6]
compression_30bar = [5]
compression_35bar = [2]
compression_38bar = [4]


dimensions = {'large': "4 1 1"}
for sys in systems:
    system_folder = f"system{sys}-8x8x25"
    system_path = base_path / system_folder
    tessilation_folder = system_path / "tessellation"
    equilibration_folder = system_path/ f"equil-{size}"

    #gets the time at avg box size
    time_at_target_box_size = avg_box_size(sys,size)

    gro_file_input = equilibration_folder / f"system{sys}-{size}-flat-extractedgro-{time_at_target_box_size}ps.gro"
    
    top_file = tessilation_folder / f"{size}-system.top"
    
    #picks the correct compression .mdp based on system composition
    if sys in compression_3bar:
        compression_mdp = mdp_path / "step6.9.1_compression.mdp"
        compression_folder = system_path /"compression"/"xzPcoupled-3bar-compression-500ns-20psreadout"
        compression_folder.mkdir(exist_ok=True)
    if sys in compression_30bar:
        compression_mdp = mdp_path / "step6.9.4_compression.mdp"
        compression_folder = system_path /"compression"/"xzPcoupled-30bar-compression-1000ns-20psreadout"
        compression_folder.mkdir(exist_ok=True)
    if sys in compression_35bar: 
        compression_mdp = mdp_path / "step6.9.2_compression.mdp"
        compression_folder = system_path /"compression"/"xzPcoupled-35bar-compression-500ns-20psreadout"
        compression_folder.mkdir(exist_ok=True)
    if sys in compression_38bar:
        compression_mdp = mdp_path / "step6.9.3_compression.mdp"
        compression_folder = system_path /"compression"/"xzPcoupled-38bar-compression-500ns-20psreadout"
        compression_folder.mkdir(exist_ok=True)

    # base filename for .tpr and .deffnm
    base_name = compression_folder / f"{size}-compression"
    compression_tpr = f"{base_name}.tpr"

    # create .tpr file
    grompp_cmd = f"gmx grompp -p {top_file} -f {compression_mdp} -c {gro_file_input} -o {compression_tpr} -maxwarn 1"
    subprocess.run(grompp_cmd, shell=True, check=True)

    # run compression
    mdrun_cmd = f"gmx mdrun -nt 24 -bonded gpu -gpu_id 1 -pin on -pinoffset 0 -v -pinstride 1 -nstlist 100 -deffnm {base_name}"
    subprocess.run(mdrun_cmd, shell=True, check=True)



 