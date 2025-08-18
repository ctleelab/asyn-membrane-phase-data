from pathlib import Path
import util as util
from util import base_path 
import subprocess

#setting paths
script_dir = Path(__file__).resolve().parent
mdp_path = Path(script_dir.parent /"mdps") 

systems = [2]
#dimensions = {'xsmall': "1 1 1", 'small': "2 1 1", 'medium': "3 1 1", 'large': "4 1 1"}
dimensions = {'large': "4 1 1"}
for sys in systems:
    system_folder = f"system{sys}-8x8x25"
    system_path = base_path / system_folder
    for size, sizes in dimensions.items():

        compression_folder = system_path / "xzPcoupled-35bar-compression-1000ns"
        compression_folder.mkdir(exist_ok=True)

        tessilation_folder = system_path / "tessilation"
        gro_file_input = tessilation_folder / f"{size}-tessilation.gro"
        # use if want to start from specific gro file
        # frame = 160140
        # gro_file_input = compression_folder / f"{size}_compression_{frame}ps.gro"
        top_file = tessilation_folder / f"{size}-system.top"
        
        #ensures you use the right compression force
        if sys == 1: 
            compression_mdp = mdp_path / "step6.8.1_compression.mdp"
        else:
            compression_mdp = mdp_path / "step6.8.2_compression.mdp"

        # base filename for .tpr and .deffnm
        base_name = compression_folder / f"{size}-compression"
        compression_tpr = f"{base_name}.tpr"

        # create .tpr file
        grompp_cmd = f"gmx grompp -p {top_file} -f {compression_mdp} -c {gro_file_input} -o {compression_tpr} -maxwarn 1"
        subprocess.run(grompp_cmd, shell=True, check=True)

        # run compression
        mdrun_cmd = f"gmx mdrun -nt 40 -update gpu -gpu_id 1 -pin off -v -pinstride 1 -nstlist 100 -deffnm {base_name}"
        subprocess.run(mdrun_cmd, shell=True, check=True)

        # #file folders for the extension run 
        # tpr_input = compression_folder/"large-compression.tpr"
        # extend_time = 1000000 

        # #extend the run
        # extend_run = f"gmx convert-tpr -s {tpr_input} -extend {extend_time} -o {compression_folder}/large-compression.tpr"
        # subprocess.run(extend_run, shell=True, check=True)

        # #extend md run 
        # extend_run_md = f"gmx mdrun -deffnm {compression_folder}/{size}-compression -nt 40 -update gpu -pin off -nstlist 100 -cpi {compression_folder}/{size}-compression.cpt"
        # subprocess.run(extend_run_md, shell=True, check=True)        




 