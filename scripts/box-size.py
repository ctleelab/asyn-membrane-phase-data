#script not set up to only do one dimension type 
#import libraries
from pathlib import Path
from util import base_path

import subprocess

import matplotlib.pyplot as plt
import numpy as np

systems = [1,2]

for sys in systems:
    #file set-up
    system_folder = f"system{sys}"
    system_path = base_path/system_folder
    equilibration_folder = f"100bar-compression"
    equilibration_path = system_path / equilibration_folder
    edr = equilibration_path / "xlarge-compression-xz-Pcoupled.edr"
    analysis_folder= system_path / f"analysis"
    analysis_folder.mkdir(exist_ok=True)
    box_size = analysis_folder / "xlarge-xzdimension-box-size-100bar-xz-Pcoupled.xvg"


    #reading edr file, to get the size of the box over time
    read_edr = f"gmx energy -f {edr} -o {box_size}"
    subprocess.run(read_edr, shell = True, check = True, input = "11\n13\n", text = True)

    xy = np.loadtxt(analysis_folder/ f"xlarge-xzdimension-box-size-100bar-xz-Pcoupled.xvg", comments = ('#', '@'))
    step = xy[:, 0]
    x_length = xy[:, 1]
    z_length = xy[:, 2]

    #create a plot
    plt.plot(step, x_length)
    plt.plot(step, z_length)
    plt.xlabel('Step')
    plt.ylabel('dimension (Nm)')
    plt.legend(['x dimension', 'z dimension'])
    plt.title(f"box x and z dimension for system{sys} at 100bar")
    plt.savefig(analysis_folder / f"xlarge-xzdimension{sys}-100bar-xz-Pcoupled.png")
    plt.close()
