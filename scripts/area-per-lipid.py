#import libraries
from pathlib import Path
from util import base_path

import subprocess

import matplotlib.pyplot as plt
import numpy as np

systems = [1,2]
lipid_number = 128

for sys in systems:
    #file set-up
    system_folder = f"system{sys}"
    system_path = base_path/system_folder
    equilibration_folder = f"equil"
    equilibration_path = system_path / equilibration_folder
    edr = equilibration_path / "equilibration6.6.edr"
    analysis_folder= system_path / f"analysis"
    analysis_folder.mkdir(exist_ok=True)
    box_size = analysis_folder / "box-size_50ns.xvg"


    #reading edr file, to get the size of the box over time
    read_edr = f"gmx energy -f {edr} -o {box_size}"
    subprocess.run(read_edr, shell = True, check = True, input = "12\n13\n\n", text = True)

    xy = np.loadtxt(analysis_folder/ f"box-size_50ns.xvg", comments = ('#', '@'))
    step = xy[:, 0]
    x_length = xy[:, 1]
    y_length = xy[:, 2]
    #calcualte the area per lipid for each step
    ApL = (x_length * y_length) / lipid_number
    #join the ApL to the original data
    np.savetxt(analysis_folder/ f'area.xvg', np.column_stack((step, ApL)))

    #create a plot
    plt.plot(step, ApL)
    plt.xlabel('Step')
    plt.ylabel('Area per lipid (nm^2)')
    plt.title(f"Area per lipid for system{sys}")

 
    plt.savefig(analysis_folder / f"area_per_lipid_system{sys}_50ns.png")
    plt.close()
