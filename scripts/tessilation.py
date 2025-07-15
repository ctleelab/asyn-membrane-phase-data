#importing libraries 
import util as util
from util import base_path 
import subprocess
import shutil

systems = [1]
dimensions = {'xsmall': "1 1 1", 'small': "2 1 1"
, 'medium': "3 1 1", 'large': "4 1 1"}


#replicating system 
for sys in systems:
    system_folder = f"system{sys}"
    system_path = base_path/system_folder
    tessilation_folder = system_path / "tessilation"
    tessilation_folder.mkdir(exist_ok =True)
    gro_file_input = system_path/"equil"/"equilibration6.7.gro"

    for size, sizes in dimensions.items():
        gro_file_output = tessilation_folder/f"{size}-tessilation.gro"
        if not gro_file_output.exists():
            tessilation_run = f"gmx genconf -f {gro_file_input} -o {gro_file_output} -nbox {sizes} -dist 0 0 0"
            subprocess.run(tessilation_run, shell=True, check=True)
        else:
            print(f"{gro_file_output} already exists.")

        #create copy of top file  
        top_file = system_path / "system.top"
        tessilation_top = tessilation_folder / f"{size}-system.top"
        if not tessilation_top.exists():
            shutil.copy(top_file, tessilation_top)
        else:
            print(f"{tessilation_top} already exists.")

        #update the top file with new tessilation number from gro file
        with open(gro_file_output, 'r') as gro_file:
            lines = gro_file.readlines()[2:]
        #set all molecule count to 0 
        DOPC_count = 0
        W = 0
        NA = 0
        CL = 0
        unique_DOPC = set()
        #go through each line in the gro file and add to the counts of each molecule type
        for line in lines:
            #only want to add if unique DOPC 1DOPC different from 2DOPC
            if "DOPC" in line:
                token = line.strip().split()[0]
                if token not in unique_DOPC:
                    unique_DOPC.add(token)
            elif "W" in line:
                W +=1
            elif "NA" in line:
                NA +=1
            elif "CL" in line:
                CL +=1
        DOPC_count = len(unique_DOPC)
        #write the updated counts to the top file
        with open(tessilation_top, 'r') as top_file:
            lines = top_file.readlines()
            unchanged_lines = lines[:9]
        with open(tessilation_top, 'w') as top_file:
            if size == 'xsmall':
                top_file.writelines(unchanged_lines)
                top_file.write(f"DOPC {DOPC_count}\n")
                top_file.write(f"W {W}\n")
                top_file.write(f"NA {NA}\n")
                top_file.write(f"CL {CL}\n")

            elif size == 'small':
                tessilations = 2
                seperated_lipid = int(DOPC_count/tessilations)
                seperated_water = int(W/tessilations)
                seperated_na = int(NA/tessilations)
                seperated_cl = int(CL/tessilations)
                top_file.writelines(unchanged_lines)
                for i in range(tessilations):
                    top_file.write(f"DOPC {seperated_lipid}\n")
                    top_file.write(f"W {seperated_water}\n")
                    top_file.write(f"NA {seperated_na}\n")
                    top_file.write(f"CL {seperated_cl}\n")
            elif size == 'medium':
                tessilations = 3
                seperated_lipid = int(DOPC_count/tessilations)
                seperated_water = int(W/tessilations)
                seperated_na = int(NA/tessilations)
                seperated_cl = int(CL/tessilations)
                top_file.writelines(unchanged_lines)
                for i in range(tessilations):
                    top_file.write(f"DOPC {seperated_lipid}\n")
                    top_file.write(f"W {seperated_water}\n")
                    top_file.write(f"NA {seperated_na}\n")
                    top_file.write(f"CL {seperated_cl}\n")
            elif size == 'large':
                tessilations = 4
                seperated_lipid = int(DOPC_count/tessilations)
                seperated_water = int(W/tessilations)
                seperated_na = int(NA/tessilations)
                seperated_cl = int(CL/tessilations)
                top_file.writelines(unchanged_lines)
                for i in range(tessilations):
                    top_file.write(f"DOPC {seperated_lipid}\n")
                    top_file.write(f"W {seperated_water}\n")
                    top_file.write(f"NA {seperated_na}\n")
                    top_file.write(f"CL {seperated_cl}\n")
            print({DOPC_count, W, NA, CL})  
