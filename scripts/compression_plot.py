#script not set up to only do one dimension type 
#import libraries
from pathlib import Path
import util
from util import base_path
from util import analysis_path
import pandas as pd

import subprocess

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import MDAnalysis as mda 
#turns the edr file into a numpy array


systems = [6]
initial_dim = "8x8x25"
final_dim = "large"
sim_time = 1000
compression = "3"
target_strain = .15

for sys in systems:
    #file set-up when files are in system folder
    system_folder = f"system{sys}-{initial_dim}"
    system_path = base_path/system_folder
    # compression_folder = f"xzPcoupled-{pressure}bar-compression-{sim_time}ns"
    #NVT compression_folder
    compression_folder = f"xzPcoupled-{compression}bar-compression-{sim_time}ns"
    compression_path = system_path / compression_folder
    edr = compression_path /"large-compression.edr"
    analysis_folder= analysis_path/f"testing_compression"
    analysis_folder.mkdir(exist_ok=True)


    aux = mda.auxiliary.EDR.EDRReader(edr)
    terms = aux.get_data(["Box-X", "Box-Z", "Pressure", "Pres-XX", "Pres-YY", "Pres-ZZ"])


    #dimension of X when time = 0
    #only do if the simulation at t=0 is flat
    time0_idx = np.argmin(np.abs(terms["Time"] - 0))
    initial_dim_x = terms["Box-X"][time0_idx]
    print(f"The initial dim: ", initial_dim_x)
    
    # # #set initial dimensionx when pulling out curvature frame
    # # initial_dim_x = 316

    compression_strain = (initial_dim_x - terms["Box-X"])/initial_dim_x
    print(compression_strain)
    Box_Xnm = terms["Box-X"]/10
    Box_Znm = terms["Box-Z"]/10
    Time_1e6 = terms["Time"]/1000000

    #combine into dataframe 
    df = pd.DataFrame({ "Time_1e6": Time_1e6,"Time": terms["Time"], "Box-X": terms["Box-X"],"Box-X-nm": Box_Xnm, "Box-Z": terms["Box-Z"],"Box-Z-nm": Box_Znm,"Compression strain": compression_strain, "Pressure": terms["Pressure"], "Pres-XX": terms["Pres-XX"], "Pres-YY": terms["Pres-YY"], "Pres-ZZ": terms["Pres-ZZ"]})

    #find the index where the compression strain is closest to the target strain
    strain_idx = np.argmin(np.abs(df["Compression strain"] - target_strain))
    time_at_target_strain = df["Time"].iloc[strain_idx]
    print(f"The strain is {target_strain} at {time_at_target_strain}ps")



    #combined plot time vs compression
    plt.plot(df["Time_1e6"], df["Compression strain"])
    plt.xlabel('time (ps)')
    plt.ylabel('compression strain')
    plt.title(f"system{sys} {final_dim}{initial_dim} at {compression}bar strain{target_strain} at {time_at_target_strain}ps ")

    # Set major ticks every 0.1 on y-axis
    plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))

    # Draw horizontal red line at compression = 0.1
    plt.axhline(y=0.1, color='red', linestyle='--', linewidth=1)
    plt.axhline(y=0.2, color='orange', linestyle='--', linewidth=1)
    plt.axhline(y=0.6, color='green', linestyle='--', linewidth=1)


    plt.tight_layout()
    plt.savefig(compression_path / f"compression_strain_system{sys}_{final_dim}{initial_dim}_{compression}bar_{sim_time}ns.png")
    plt.close()
