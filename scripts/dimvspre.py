#Goal plot the x-dimension (length of the box in the x-dimension) vs pressure over the simulation 
#for now I will use the 6x6x12nm (large) system with DOPC or DPPC as the membrane appears to have one local maximum 


#for the x-dimension and x-pressure I would need the .edr file

#import libraries 
from pathlib import Path
from util import base_path

import matplotlib.pyplot as plt
import numpy as np

import MDAnalysis as mda

initial_dim = "6x6x12"
#dimension and total lipid number 
final_dims = ['large']

#will be created in current working directory, which should be scipts
figures = Path("figures")
figures.mkdir(parents=True, exist_ok=True)
figures_dimvspre = figures/"dimvspre"
figures_dimvspre.mkdir(parents=True, exist_ok=True)

# Define systems and corresponding line styles and labels
systems = [1, 2]
line_styles = {1: '-', 2: '--'}
labels = {1: "DOPC", 2: "DPPC"}

for dim in final_dims:
    # Create one figure for this dim, comparing system1 and system2
    fig, axs = plt.subplots(3, figsize=(7, 8))

    for sys in systems:
        # File setup
        system_folder = f"system{sys}-{initial_dim}"
        system_path = base_path / system_folder
        if sys == 1: 
            pressure_folder = "100bar-xzPcoupled-compression"
        else:
            pressure_folder = "200bar-xzPcoupled-compression"

        pressure_path = system_path / pressure_folder

        file_name = f"{dim}-compression"
        edr = pressure_path / f"{file_name}.edr"

        # Read EDR data
        aux = mda.auxiliary.EDR.EDRReader(edr)
        terms = aux.get_data(["Box-X", "Box-Z", "Pres-XX", "Time"])

        # Plot system data with corresponding line style
        axs[0].plot(terms["Time"], terms["Box-X"], linestyle=line_styles[sys], label=labels[sys])
        axs[1].plot(terms["Time"], terms["Box-Z"], linestyle=line_styles[sys], label=labels[sys])
        axs[2].plot(terms["Time"], terms["Pres-XX"], linestyle=line_styles[sys], label=labels[sys])

    # Set titles and labels
    axs[0].set_title("Box-X dimension")
    axs[1].set_title("Box-Z dimension")
    axs[2].set_title("Pressure-XX")

    for ax in axs:
        ax.set_xlabel("Time (ps)")
        ax.legend()

    plt.tight_layout()
    plt.savefig(figures_dimvspre / f"compare-system1-2-{initial_dim}-{dim}-100vs200bar.png")
    plt.close()
