#script not set up to only do one dimension type 
#import libraries
from pathlib import Path
from util import base_path
import pandas as pd

import subprocess

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda 
#turns the edr file into a numpy array


systems = [5]
size = "medium"

for sys in systems:
    #file set-up
    system_folder = f"system{sys}-8x8x25"
    system_path = base_path/system_folder
    #compression_folder = f"xzPcoupled-compression"
    equilibration_folder = f"equil-{size}"
    equilibration_path = system_path / equilibration_folder
    edr = equilibration_path / "equilibration6.7.edr"


    aux = mda.auxiliary.EDR.EDRReader(edr)
    terms = aux.get_data(["Box-X", "Box-Z", "Time"])

    df = pd.DataFrame({"Time": terms["Time"], "Box-X": terms["Box-X"], "Box-Z": terms["Box-Z"]})
    box_z_avg = np.mean(terms["Box-Z"])
    box_x_avg = np.mean(terms["Box-X"])
    

    #determine the frame where the intial_dim_x is the average 
    target_xdim_idx = np.argmin(np.abs(df["Box-X"] - box_x_avg))

    time_at_target_box_size = df["Time"].iloc[target_xdim_idx]


    #create a plot
    plt.plot(terms["Time"],terms["Box-X"])
    plt.xlabel('Time')
    plt.ylabel('dimension (Nm)')
    plt.legend('x dimension')
    plt.title(f"box x dimension average-Xdim {box_x_avg} at time {time_at_target_box_size}")
    plt.axhline(box_x_avg, linestyle = "--", linewidth = 2)
    plt.savefig(equilibration_path / f"{size}-x-box-dimension-sys{sys}-equil6.7.png")
    plt.close()


    #create a plot
    plt.plot(terms["Time"],terms["Box-Z"])
    plt.xlabel('Time')
    plt.ylabel('dimension (Nm)')
    plt.legend('z dimension')
    plt.title(f"box  z dimension for system{sys} average-Zdim {box_z_avg}")
    plt.axhline(box_z_avg, linestyle = "--", linewidth = 2)
    plt.savefig(equilibration_path / f"{size}-z-box-dimension-sys{sys}-equil6.7.png")
    plt.close()