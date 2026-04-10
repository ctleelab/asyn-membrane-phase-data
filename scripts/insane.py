#helps create the initial system (8x8x25nm) that later gets tessellated to generate the larger systems. 

import util as util
from util import base_path
from pathlib import Path
import subprocess 

#defining script directory
script_dir = Path(__file__).resolve().parent
martini_itp_dir = script_dir.parent / "martini3-files"

#stating completed systems
systems = [1,2,3,4,5,6]

for sims, compositions in util.system_compositions.items():
    if sims not in systems:
       continue 
    
    #create a folder for each system in the systems folder
    folder_name = f"system{sims}-8x8x25/min-xsmall"
    #coverts the folder name which is a string to a Path object 
    folder_path = base_path/folder_name
    #creates the folder 
    folder_path.mkdir(exist_ok=True)
    
    #turn the item in the dictionary (compositions) into a string for insane to properly run
    for k, v in compositions.items():
        compositions_str = f"{k}:{v}"
    
    #write gro and top files to the created folder with no neutralization
    output_gro = folder_path / f"system-{sims}-8nm.gro"
    output_top = folder_path / f"system.top"
    output_tpr = folder_path / f"system-{sims}-8nm.tpr"


    #write gro and top files to the created folder 
    #this naming need to match the CS_minimization.sh script
    neutral_gro = folder_path / f"neutral.gro"

    #run insane function 
    insane_run = f"insane -x 8 -y 8 -z 25 -sol W -o {output_gro} -l {compositions_str} -p {output_top}"

    #subprocess running, meaning takes the insane_run string and runs it as a bash command
    subprocess.run(insane_run ,shell = True, check = True )

    #rewrite the phrase #include "martini.itp" with #include "/scratch/casakurai/asyn-phase-binding-data/toppar/martini_v3.0_openbeta/martini.itp"
    with open(output_top, 'r') as file:
        lines = file.readlines()
    new_lines = []
    for line in lines:
        if line.strip().startswith('#include "martini.itp"'):
            new_lines.append(f'#include "{martini_itp_dir}/martini.itp"\n')
        else:
            new_lines.append(line)  # Keep all other lines unchanged


    # Write back the modified .top file
    with open(output_top, 'w') as f:
        f.writelines(new_lines)


    #create a tpr file with step6.0_minimization.mdp in order to neutralize the membrane 
    initial_run = f"gmx grompp -f ../mdps/step6.0_minimization.mdp -c {output_gro} -p {output_top} -o {output_tpr}" 

    
    subprocess.run(initial_run, shell=True, check=True)

    #neutralize the membrane with .15 M NaCl
    neutralize_run = f"gmx genion -s {output_tpr} -o {neutral_gro} -p {output_top} -neutral -conc .15"

    #input ensures that Waters is replaced by Na and Cl ions
    subprocess.run(neutralize_run, shell = True, check = True, input = 'W\n', text = True)


    