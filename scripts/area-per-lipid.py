#import libraries
from pathlib import Path
from util import base_path
from util import analysis_path
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda

systems = [1]
lipid_number = 800
initial_dim = "8x8x25"
final_dim = "large"
pressure_initial = "2.25"
sim_time = "400"

for sys in systems:
    #file set-up
    system_folder = f"system{sys}-{initial_dim}"
    system_path = base_path/system_folder
    pressure_folder = f"xzPcoupled-{pressure_initial}bar-compression-{sim_time}ns"
    pressure_path = system_path/pressure_folder
    edr = pressure_path/f"{final_dim}-compression.edr"

    #code for edr file in equilibration folder
    # equilibration_folder = f"equil"
    # equilibration_path = system_path / equilibration_folder
    # edr = equilibration_path / "equilibration6.6.edr"


    analysis_folder= analysis_path/f"testing_compression"
    analysis_folder.mkdir(exist_ok=True)

    aux = mda.auxiliary.EDR.EDRReader(edr)
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

 
    plt.savefig(analysis_folder / f"area_per_lipid_system{sys}_{final_dim}{initial_dim}_{pressure_initial}bar_{sim_time}ns.png")
    plt.close()
