################
##GENERAL INFO##
################
#can only run this script for either system 1 or 2. Cannot do at the same time as it counts either DOPC or DPPC 

#importing libraries 
import util as util
from util import base_path 
import subprocess
import shutil

systems = [3]
#dimensions = {'xsmall': "1 1 1", 'small': "2 1 1", 'medium': "3 1 1", 'large': "4 1 1", 'xlarge': "10 1 1"}
dimensions = {'large': "4 1 1"}
lipids = ["DOPC", "DOPS", "DOPA", "DPPC", "DPPS", "DPPA"]
ions = ["W", "CL", "NA"]

#replicating system 
for sys in systems:
    system_folder = f"system{sys}-8x8x25"
    system_path = base_path/system_folder
    tessilation_folder = system_path / "tessellation"
    tessilation_folder.mkdir(exist_ok =True)
    gro_file_input = system_path/"equil-xsmall"/"equilibration6.7.gro"

    for size, sizes in dimensions.items():
        gro_file_output = tessilation_folder/f"{size}-tessilation.gro"
        if not gro_file_output.exists():
            tessilation_run = f"gmx genconf -f {gro_file_input} -o {gro_file_output} -nbox {sizes} -dist 0 0 0"
            subprocess.run(tessilation_run, shell=True, check=True)
        else:
            print(f"{gro_file_output} already exists.")

        #create copy of top file  
        top_file = system_path /"min-xsmall" /"system.top"
        tessilation_top = tessilation_folder / f"{size}-system.top"
        if not tessilation_top.exists():
            shutil.copy(top_file, tessilation_top)
        else:
            print(f"{tessilation_top} already exists.")

        #update the top file with new tessilation number from gro file
        with open(gro_file_output, 'r') as gro_file:
            lines = gro_file.readlines()[2:]
        
        #set all molecules counts as 0
        W = 0
        NA = 0
        CL = 0
        unique_DOPC = set()
        unique_DPPC = set()
        unique_DOPA = set()
        unique_DPPA = set()
        unique_DOPS = set()
        unique_DPPS = set()
        for line in lines:
            if "DOPC" in line:
                token = line.strip().split()[0]
                if token not in unique_DOPC:
                    unique_DOPC.add(token)
            if "DPPC" in line:
                token = line.strip().split()[0]
                if token not in unique_DPPC:
                    unique_DPPC.add(token)
            if "DOPA" in line:
                token = line.strip().split()[0]
                if token not in unique_DOPA:
                    unique_DOPA.add(token)
            if "DPPA" in line:
                token = line.strip().split()[0]
                if token not in unique_DPPA:
                    unique_DPPA.add(token)
            if "DOPS" in line:
                token = line.strip().split()[0]
                if token not in unique_DOPS:
                    unique_DOPS.add(token)
            if "DPPS" in line:
                token = line.strip().split()[0]
                if token not in unique_DPPS:
                    unique_DPPS.add(token)
            if "W" in line:
                W +=1
            if "NA" in line:
                NA +=1
            if "CL" in line:
                CL +=1 
        DOPC_count = len(unique_DOPC)
        DPPC_count = len(unique_DPPC)
        DOPA_count = len(unique_DOPA)
        DPPA_count = len(unique_DPPA)
        DOPS_count = len(unique_DOPS)
        DPPS_count = len(unique_DPPS)

        # #covert to counts per ion or lipid type
        # lipid_counts = {lip: len(resids) for lip, resids in lipid_residues.items()}
        
        # W = ion_counts["W"]
        # NA = ion_counts["NA"]
        # CL = ion_counts["CL"]
                    

        #write the updated counts to the .top file
        #but write in a format that matches the .gro file
        with open(tessilation_top, 'r') as top_file:
            lines = top_file.readlines()
            unchanged_lines = lines[:9]
        with open(tessilation_top, 'w') as top_file:
            if size == 'xsmall':
                top_file.writelines(unchanged_lines)
                tessellations = 1
                seperated_DOPC = DOPC_count/tessellations
                seperated_DOPS = DOPS_count/tessellations
                seperated_DOPA = DOPA_count/tessellations
                seperated_DPPC = DPPC_count/tessellations
                seperated_DPPS = DPPS_count/tessellations
                seperated_DPPA = DPPA_count/tessellations
                seperated_W = W/tessellations
                seperated_NA = NA/tessellations
                seperated_CL = CL/tessellations
                # lipids (only write if present)
                if DOPC_count > 0:
                    top_file.write(f"DOPC {seperated_DOPC}\n")
                if DOPS_count > 0:
                    top_file.write(f"DOPS {seperated_DOPS}\n")
                if DOPA_count > 0:
                    top_file.write(f"DOPA {seperated_DOPA}\n")
                if DPPC_count > 0:
                    top_file.write(f"DPPC {seperated_DPPC}\n")
                if DPPS_count > 0:
                    top_file.write(f"DPPS {seperated_DPPS}\n")
                if DPPA_count > 0:
                    top_file.write(f"DPPA {seperated_DPPA}\n")
                top_file.write(f"W {seperated_W}\n")
                top_file.write(f"NA {seperated_NA}\n")
                top_file.write(f"CL {seperated_CL}\n")


            elif size == 'small':
                tessellations = 2
                
                seperated_DOPC = (DOPC_count/tessellations)/2 #divide by two because two leaflets
                seperated_DOPS = (DOPS_count/tessellations)/2
                seperated_DOPA = (DOPA_count/tessellations)/2
                seperated_DPPC = (DPPC_count/tessellations)/2
                seperated_DPPS = (DPPS_count/tessellations)/2
                seperated_DPPA = (DPPA_count/tessellations)/2
                seperated_W = W/tessellations
                seperated_NA = NA/tessellations
                seperated_CL = CL/tessellations
                top_file.writelines(unchanged_lines)
                for i in range(tessellations):
                    if DOPC_count > 0:
                        top_file.write(f"DOPC {seperated_DOPC}\n")
                    if DOPS_count > 0:
                        top_file.write(f"DOPS {seperated_DOPS}\n")
                    if DOPA_count > 0:
                        top_file.write(f"DOPA {seperated_DOPA}\n")
                    if DPPC_count > 0:
                        top_file.write(f"DPPC {seperated_DPPC}\n")
                    if DPPS_count > 0:
                        top_file.write(f"DPPS {seperated_DPPS}\n")
                    if DPPA_count > 0:
                        top_file.write(f"DPPA {seperated_DPPA}\n")

                    top_file.write(f"W {seperated_W}\n")
                    top_file.write(f"NA {seperated_NA}\n")
                    top_file.write(f"CL {seperated_CL}\n")

            elif size == 'medium':
                tessellations = 3
                seperated_DOPC = (DOPC_count/tessellations)/2 #divide by two because two leaflets
                seperated_DOPS = (DOPS_count/tessellations)/2
                seperated_DOPA = (DOPA_count/tessellations)/2
                seperated_DPPC = (DPPC_count/tessellations)/2
                seperated_DPPS = (DPPS_count/tessellations)/2
                seperated_DPPA = (DPPA_count/tessellations)/2
                seperated_W = W/tessellations
                seperated_NA = NA/tessellations
                seperated_CL = CL/tessellations
                top_file.writelines(unchanged_lines)
                for i in range(tessellations):
                    if DOPC_count > 0:
                        top_file.write(f"DOPC {seperated_DOPC}\n")
                        print(seperated_DOPC)
                    if DOPS_count > 0:
                        top_file.write(f"DOPS {seperated_DOPS}\n")
                        print(seperated_DOPS)
                    if DOPA_count > 0:
                        top_file.write(f"DOPA {seperated_DOPA}\n")
                    if DPPC_count > 0:
                        top_file.write(f"DPPC {seperated_DPPC}\n")
                    if DPPS_count > 0:
                        top_file.write(f"DPPS {seperated_DPPS}\n")
                    if DPPA_count > 0:
                        top_file.write(f"DPPA {seperated_DPPA}\n")

                    top_file.write(f"W {seperated_W}\n")
                    top_file.write(f"NA {seperated_NA}\n")
                    top_file.write(f"CL {seperated_CL}\n")

            elif size == 'large':
                tessellations = 4
                seperated_DOPC = (DOPC_count/tessellations)/2 #divide by two because two leaflets
                seperated_DOPS = (DOPS_count/tessellations)/2
                seperated_DOPA = (DOPA_count/tessellations)/2
                seperated_DPPC = (DPPC_count/tessellations)/2
                seperated_DPPS = (DPPS_count/tessellations)/2
                seperated_DPPA = (DPPA_count/tessellations)/2
                seperated_W = W/tessellations
                seperated_NA = NA/tessellations
                seperated_CL = CL/tessellations
                top_file.writelines(unchanged_lines)
                for i in range(tessellations):
                    if DOPC_count > 0:
                        top_file.write(f"DOPC {seperated_DOPC}\n")
                    if DOPS_count > 0:
                        top_file.write(f"DOPS {seperated_DOPS}\n")
                    if DOPA_count > 0:
                        top_file.write(f"DOPA {seperated_DOPA}\n")
                    if DPPC_count > 0:
                        top_file.write(f"DPPC {seperated_DPPC}\n")
                    if DPPS_count > 0:
                        top_file.write(f"DPPS {seperated_DPPS}\n")
                    if DPPA_count > 0:
                        top_file.write(f"DPPA {seperated_DPPA}\n")

                    top_file.write(f"W {seperated_W}\n")
                    top_file.write(f"NA {seperated_NA}\n")
                    top_file.write(f"CL {seperated_CL}\n")
            elif size == 'xlarge-z' :
                tessellations = 50
                seperated_DOPC = (DOPC_count/tessellations)/2 #divide by two because two leaflets
                seperated_DOPS = (DOPS_count/tessellations)/2
                seperated_DOPA = (DOPA_count/tessellations)/2
                seperated_DPPC = (DPPC_count/tessellations)/2
                seperated_DPPS = (DPPS_count/tessellations)/2
                seperated_DPPA = (DPPA_count/tessellations)/2
                seperated_W = W/tessellations
                seperated_NA = NA/tessellations
                seperated_CL = CL/tessellations
                top_file.writelines(unchanged_lines)
                for i in range(tessellations):
                    if DOPC_count > 0:
                        top_file.write(f"DOPC {seperated_DOPC}\n")
                    if DOPS_count > 0:
                        top_file.write(f"DOPS {seperated_DOPS}\n")
                    if DOPA_count > 0:
                        top_file.write(f"DOPA {seperated_DOPA}\n")
                    if DPPC_count > 0:
                        top_file.write(f"DPPC {seperated_DPPC}\n")
                    if DPPS_count > 0:
                        top_file.write(f"DPPS {seperated_DPPS}\n")
                    if DPPA_count > 0:
                        top_file.write(f"DPPA {seperated_DPPA}\n")

                    top_file.write(f"W {seperated_W}\n")
                    top_file.write(f"NA {seperated_NA}\n")
                    top_file.write(f"CL {seperated_CL}\n")
            elif size == 'xlarge' :
                tessellations = 10
                seperated_DOPC = (DOPC_count/tessellations)/2 #divide by two because two leaflets
                seperated_DOPS = (DOPS_count/tessellations)/2
                seperated_DOPA = (DOPA_count/tessellations)/2
                seperated_DPPC = (DPPC_count/tessellations)/2
                seperated_DPPS = (DPPS_count/tessellations)/2
                seperated_DPPA = (DPPA_count/tessellations)/2
                seperated_W = W/tessellations
                seperated_NA = NA/tessellations
                seperated_CL = CL/tessellations
                top_file.writelines(unchanged_lines)
                for i in range(tessellations):
                    if DOPC_count > 0:
                        top_file.write(f"DOPC {seperated_DOPC}\n")
                    if DOPS_count > 0:
                        top_file.write(f"DOPS {seperated_DOPS}\n")
                    if DOPA_count > 0:
                        top_file.write(f"DOPA {seperated_DOPA}\n")
                    if DPPC_count > 0:
                        top_file.write(f"DPPC {seperated_DPPC}\n")
                    if DPPS_count > 0:
                        top_file.write(f"DPPS {seperated_DPPS}\n")
                    if DPPA_count > 0:
                        top_file.write(f"DPPA {seperated_DPPA}\n")

                    top_file.write(f"W {seperated_W}\n")
                    top_file.write(f"NA {seperated_NA}\n")
                    top_file.write(f"CL {seperated_CL}\n")
