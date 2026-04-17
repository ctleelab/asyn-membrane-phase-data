#import libraries
from pathlib import Path
from util import base_path
import pandas as pd

import subprocess

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda 
#turns the edr file into a numpy array


def avg_box_size(sys,size):
    '''
    Determines the average box size in the x dimension. 
    Plots the box size over time.
    Then extracts the frame closest to the average box size and saves it as a .gro file. 

    Args
    sys: int, the system number to analyze
    size: str, the size of the system being analyzed.

    Returns
    Extract gro frame closest to avg Lx box size
    Plot of box size over time in Lx and Lz

    '''

    #file set-up
    system_folder = f"system{sys}-8x8x25"
    system_path = base_path/system_folder
    #compression_folder = f"xzPcoupled-compression"
    equilibration_folder = f"equil-{size}"
    equilibration_path = system_path / equilibration_folder
    edr = equilibration_path / "equilibration6.9-ext.edr"


    aux = mda.auxiliary.EDR.EDRReader(edr)
    terms = aux.get_data(["Box-X", "Box-Z", "Time"])

    df = pd.DataFrame({"Time": terms["Time"], "Box-X": terms["Box-X"], "Box-Z": terms["Box-Z"]})
    box_z_avg = np.mean(terms["Box-Z"])
    box_x_avg = np.mean(terms["Box-X"])
    

    #determine the frame where the intial_dim_x is the average 
    target_xdim_idx = np.argmin(np.abs(df["Box-X"] - box_x_avg))
    #determine the time of the target frame 
    time_at_target_box_size = df["Time"].iloc[target_xdim_idx]
    actual_xdim = df["Box-X"].iloc[target_xdim_idx]


    #creates a plot of Lx over time 
    plt.plot(terms["Time"],terms["Box-X"])
    plt.xlabel('Time')
    plt.ylabel('dimension (A)')
    plt.legend('x dimension')
    plt.title(f"box x dimension average-Xdim {box_x_avg:.5f}. At time {time_at_target_box_size} is {actual_xdim:.5f} ")
    plt.axhline(box_x_avg, linestyle = "--", linewidth = 2)
    plt.savefig(equilibration_path / f"{size}-x-box-dimension-sys{sys}-equil6.9-ext.png")
    plt.close()

    #creates a plot of Lx over time 
    plt.plot(terms["Time"],terms["Box-Z"])
    plt.xlabel('Time')
    plt.ylabel('dimension (A)')
    plt.legend('z dimension')
    plt.title(f"box  z dimension for system{sys} average-Zdim {box_z_avg}.")
    plt.axhline(box_z_avg, linestyle = "--", linewidth = 2)
    plt.savefig(equilibration_path / f"{size}-z-box-dimension-sys{sys}-equil6.9-ext.png")
    plt.close()


    #extract gro frame closest to avg Lx box size
    output_file = f"system{sys}-{size}-avg-flat-extractedgro"
    file_name = f"equilibration6.9-ext"
    output_path = Path(f"{equilibration_path}/{output_file}.gro")
    if not output_path.exists():
        frame_trr = (
            f"echo '0' | gmx trjconv "
            f"-f {equilibration_path}/{file_name}.xtc "
            f"-s {equilibration_path}/{file_name}.tpr "
            f"-o {output_path} "
            f"-dump {time_at_target_box_size}"
        )
        subprocess.run(frame_trr, shell=True, check=True)
    else:
        print(f"File already exists: {output_path}")
    frame_trr = f"echo '0'| gmx trjconv -f {equilibration_path}/{file_name}.xtc -s {equilibration_path}/{file_name}.tpr -o {equilibration_path}/{output_file}.gro -dump {time_at_target_box_size}"

    return time_at_target_box_size

