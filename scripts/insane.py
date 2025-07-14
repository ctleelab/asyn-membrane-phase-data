#load in the systems 
import util as util
from util import base_path

#create a folder
from pathlib import Path

#helps run bash commands within python 
import subprocess 


#goal loop through each system composition in util 
#create a for loop takes the component system compositions and iterates through the dictionary 

#defining script directory
#_file_ is the path where the python script is currently being exectuted
#.resolve gives you the absolute path 
#.parent takes you to the directory where the script is contained 
script_dir = Path(__file__).resolve().parent
#goes up one folder from the parent, and then into martini3-files
martini_itp_dir = script_dir.parent / "martini3-files"

# #stating completed systems
completed_systems = [1,2]

#items because system compositions are a dictionary
for sims, compositions in util.system_compositions.items():
    if sims not in completed_systems:
       continue 
    
    #create a folder for each system in the systems folder
    folder_name = f"system{sims}"
    #coverts the folder name which is a string to a Path object 
    folder_path = base_path/folder_name
    #creates the folder 
    folder_path.mkdir(exist_ok=True)
    
    #turn the item in the dictionary (compositions) into a string for insane to properly run
    for k, v in compositions.items():
        compositions_str = f"{k}:{v}"
    
    #write gro and top files to the created folder with no neutralization
    output_gro = folder_path / f"system-{sims}-6nm.gro"
    output_top = folder_path / f"system.top"
    output_tpr = folder_path / f"system-{sims}-6nm.tpr"


    #write gro and top files to the created folder 
    #this naming need to match the CS_minimization.sh script
    neutral_gro = folder_path / f"neutral.gro"

    #run insane function for large box 
    insane_run = f"insane -x 6 -y 6 -z 6 -sol W -o {output_gro} -l {compositions_str} -p {output_top}"

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

    #neutralize the membrane with NaCl
    neutralize_run = f"gmx genion -s {output_tpr} -o {neutral_gro} -p {output_top} -neutral -conc .15"

    #input ensures that Waters is replaced by Na and Cl ions
    subprocess.run(neutralize_run, shell = True, check = True, input = 'W\n', text = True)


    