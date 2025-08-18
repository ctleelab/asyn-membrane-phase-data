#note this won't work on the NVT ensemble, because the volume does not fluctuate 
#import libraries
from pathlib import Path
from util import base_path
from util import analysis_path
import subprocess
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda

#parameters
systems = [1]
lipid_number = 800
initial_dim = "8x8x25"
final_dim = "large"
pressure_initial = 3
sim_time = 200
frame = 85620
equil_step = 7.2
lipid = "DPPC"

for sys in systems:
    #file set-up
    analysis_folder= analysis_path/f"curvature_selection"
    pressure_folder = analysis_folder/f"{pressure_initial}bar_{sim_time}ns"

    system_folder = f"system{sys}-{initial_dim}-{lipid}-{frame}ps-NVT"
    curvature_folder = f"system{sys}-{initial_dim}-{lipid}-{frame}ps-NVT"
    system_path = pressure_folder/system_folder/curvature_folder
    edr = system_path/"equil"/f"equilibration{equil_step}.edr"


    #where files will save 
    ApL_analysis = pressure_folder/f"ApL_plot"
    ApL_analysis.mkdir(exist_ok = True)

    aux = mda.auxiliary.EDR.EDRReader(edr)

    print(aux.terms)
    terms = aux.get_data(["Box-X", "Box-Y", "Time"])
    box_x = aux.get_data("Box-X")
    box_y = aux.get_data("Box-Y")
    time = box_x["Time"]

    #calculate the area per lipid for each step
    ApL = (np.array(box_x["Box-X"]) * np.array(box_y["Box-Y"])) / lipid_number
    
    #join the ApL to the original data

    #create a plot
    plt.plot(time, ApL)
    plt.xlabel('time(ps)')
    plt.ylabel('Area per lipid (nm^2)')
    plt.title(f"Area per lipid for system{sys}")

 
    plt.savefig(ApL_analysis/f"area_per_lipid_system{sys}_{final_dim}{initial_dim}_{pressure_initial}bar_{sim_time}ns.png")
    plt.close()



    #dimension of X when time = 0
    #the Lo needs to be the starting length of minimization box
    time0_idx = np.argmin(np.abs(terms["Time"] - 0))
    initial_dim_x = terms["Box-X"][time0_idx]
    print(initial_dim_x)

    compression_strain = (initial_dim_x- terms["Box-X"])/initial_dim_x
    Box_Xnm = terms["Box-X"]/10

    #combine into dataframe 
    df = pd.DataFrame({"Time": terms["Time"], "Box-X": terms["Box-X"],"Box-X-nm": Box_Xnm, "Box-Z": terms["Box-Z"],"Compression strain": compression_strain})

    #combined plot time vs compression
    plt.plot(df["Time"], df["Compression strain"])
    plt.xlabel('time(ps)')
    plt.ylabel('compression strain')
    plt.title(f"compression strain for system{sys} {final_dim}{initial_dim} at {pressure_initial}bar for {sim_time}ns")
    plt.savefig(ApL_analysis / f"compression_strain_system{sys}_{final_dim}{initial_dim}_{pressure_initial}bar_{sim_time}ns_frame{frame}_{compression_pressure}barx.png")
    plt.close()
