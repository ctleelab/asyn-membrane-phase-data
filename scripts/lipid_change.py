#change DOPC to DPPC .gro
#change DOPC to DPPC
#change D2A to C2A
#change D2B to C2B 

import util
from util import base_path 
from util import analysis_path
from pathlib import Path 
import subprocess
import MDAnalysis as mda
import shutil

#setting paths
sys = 1
#use if you are looking directly at the vmd simulation 
# vmd_frame = 50
# frame =  int(vmd_frame*1000*.02)
#use if you're pulling from a particular time (ps)
frame =  int(85620)
initial_dim = "8x8x25"
final_dim = "large"
compression = "3bar"
lipid = "DPPC"
sim_time = 200
condition = "NVT"

#input files 
system_folder = f"system{sys}-{initial_dim}"
compresssion_folder = f"xzPcoupled-{compression}-compression-{sim_time}ns-original"
input_folder = base_path/system_folder/compresssion_folder
gro_file = input_folder/f"{final_dim}_compression_{frame}ps.gro"
system_top_folder = f"system{sys}-{initial_dim}"
top_file = base_path/system_top_folder/"tessilation"/f"{final_dim}-system.top"

#where files are saving
system_analysis_time_folder = f"system{sys}-{initial_dim}-{lipid}-{frame}ps"
analysis = util.analysis_path
analysis_curv_folder = f"curvature_selection"
system_analysis_folder = analysis/analysis_curv_folder/f"{compression}_{sim_time}ns"/f"system{sys}-{initial_dim}-{lipid}-{frame}ps-{condition}"
system_analysis_folder.mkdir(exist_ok = True)
analysis_curv_dir = analysis/analysis_curv_folder/f"{compression}"/f"{system_analysis_folder}"/f"{system_analysis_time_folder}"
analysis_curv_dir.mkdir(exist_ok = True)


#goal: change DOPC to DPPC in top files
#first copy file 
top_file_copy = analysis_curv_dir/f"{final_dim}-system.top"

with open(top_file, "r") as f:
    content = f.read()

content = content.replace("DOPC", lipid)
with open(top_file_copy, "w") as f:
    f.write(content)



#goal: change DOPC to DPPC in .gro files
#first copy file 
if lipid == "DPPC":
    gro_file_copy = analysis_curv_dir/f"{final_dim}_compression_{frame}ps.gro"
    with open(gro_file, "r") as f:
        content = f.read()

    content = content.replace("DOPC", lipid).replace("D2A", "C2A").replace("D2B", "C2B")

    with open(gro_file_copy, "w") as f:
        f.write(content)
else:
    gro_file_copy = analysis_curv_dir/f"{final_dim}_compression_{frame}ps.gro"
    shutil.copy(gro_file, gro_file_copy)

