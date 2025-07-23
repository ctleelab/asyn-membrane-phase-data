from pathlib import Path
import util as util
from util import base_path 
import subprocess

#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps") 

systems = [2]
#dimensions = {'xsmall': "1 1 1", 'small': "2 1 1", 'medium': "3 1 1", 'large': "4 1 1"}
dimensions = {'medium': "3 1 1"}
for sys in systems:
    system_folder = f"system{sys}-12.5x12.5x25"
    system_path = base_path / system_folder
    for size, sizes in dimensions.items():

        compression_folder = system_path / "xzPcoupled-100bar-compression"
        compression_folder.mkdir(exist_ok=True)

        tessilation_folder = system_path / "tessilation"
        gro_file_input = tessilation_folder / f"{size}-tessilation.gro"
        top_file = tessilation_folder / f"{size}-system.top"
        compression_mdp = mdp_path / "step6.8_compression.mdp"

        # base filename for .tpr and .deffnm
        base_name = compression_folder / f"{size}-compression"
        compression_tpr = f"{base_name}.tpr"

        # create .tpr file
        grompp_cmd = f"gmx grompp -p {top_file} -f {compression_mdp} -c {gro_file_input} -o {compression_tpr} -maxwarn 1"
        subprocess.run(grompp_cmd, shell=True, check=True)

        # run compression
        mdrun_cmd = f"gmx mdrun -nt 32 -update gpu -gpu_id 0 -pin on -v -pinstride 1 -nstlist 100 -deffnm {base_name}"
        subprocess.run(mdrun_cmd, shell=True, check=True)
        




 