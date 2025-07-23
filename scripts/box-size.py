#script not set up to only do one dimension type 
#import libraries
from pathlib import Path
from util import base_path

import subprocess

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda 
#turns the edr file into a numpy array


systems = [1,2]

for sys in systems:
    #file set-up
    system_folder = f"system{sys}-12.5x12.5x25"
    system_path = base_path/system_folder
    compression_folder = f"xzPcoupled-compression"
    equilibration_path = system_path / compression_folder
    edr = equilibration_path / "large-compression.edr"
    analysis_folder= system_path / f"analysis"
    analysis_folder.mkdir(exist_ok=True)


    aux = mda.auxiliary.EDR.EDRReader(edr)
    terms = aux.get_data(["Box-X", "Box-Z"])


    #create a plot
    plt.plot(terms["Time"],terms["Box-X"])
    plt.plot(terms["Time"],terms["Box-Z"])
    plt.xlabel('Time')
    plt.ylabel('dimension (Nm)')
    plt.legend(['x dimension', 'z dimension'])
    plt.title(f"box x and z dimension for system{sys} at 200bar")
    plt.savefig(analysis_folder / f"large-xzdimension{sys}-200bar-xz-Pcoupled.png")
    plt.close()
