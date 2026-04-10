# import libraries
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import MDAnalysis as mda

# turns the edr file into a numpy array


systems = [1,2,3,4,5,6]
target_strain = 0.2
# compare strain using Lo from flat membranes vs Lo of buckled
system_strain_comparison = {}
compression_3bar = [1,3,6]
compression_30bar = [5]
compression_35bar = [2]
compression_38bar = [4]


def compression_plot_set_Lo(system):
    '''
    Determine the time at which the compression strain is closest to the target strain.

    Args
    system: int, the system number to analyze

    Returns
    df: pandas dataframe with columns Time_1e6, Time, Box-X, Compression strain
    compression_at_target_strain: float, the actual compression strain at the target strain
    time_at_target_strain: float, the time at which the the strain is the closest to the target strain
    specific_compression_folder: Path, the path to the specific compression folder for the system

    '''
    # first grab the length of Lo which is the start of the compression
    compression_folder = Path(
        f"/scratch/local/casakurai/asyn-phase-binding-data/systems/system{system}-8x8x25/compression/"
    )
    if system in compression_3bar:
        specific_compression_folder = compression_folder/"xzPcoupled-3bar-compression-500ns-20psreadout"
    if system in compression_30bar: 
        specific_compression_folder = compression_folder/"xzPcoupled-30bar-compression-1000ns-20psreadout"
    if system in compression_35bar: 
        specific_compression_folder = compression_folder/"xzPcoupled-35bar-compression-500ns-20psreadout"
    if system in compression_38bar:
        specific_compression_folder = compression_folder/"xzPcoupled-38bar-compression-500ns-20psreadout"



    edr_file = specific_compression_folder/"large-compression.edr"
    aux_flat = mda.auxiliary.EDR.EDRReader(edr_file)

    terms_compression = aux_flat.get_data(["Box-X", "Box-Z", "Time"])


    #get box-x at time = 0

    time0_idx = np.argmin(np.abs(terms_compression["Time"] - 0))
    Lo  = terms_compression["Box-X"][time0_idx]


    

    # from buckled
    compression_strain = (
        Lo - terms_compression["Box-X"]
    ) / Lo
    Time_1e6 = terms_compression["Time"] / 1000000

    # combine into dataframe
    df = pd.DataFrame(
        {
            "Time_1e6": Time_1e6,
            "Time": terms_compression["Time"],
            "Box-X": terms_compression["Box-X"],
            "Compression strain": compression_strain,
        }
    )

    # find the index where the compression strain is closest to the target strain
    # this might not be entirely accurate
    strain_idx = np.argmin(np.abs(df["Compression strain"] - target_strain))

    time_at_target_strain = df["Time"].iloc[strain_idx]
    compression_at_target_strain = df["Compression strain"].iloc[strain_idx]
    Lt_at_target_strain = df["Box-X"].iloc[strain_idx]


    system_strain_comparison[system] = [
        target_strain,
        compression_at_target_strain
    ]

    return df,compression_at_target_strain, time_at_target_strain,specific_compression_folder


def compression_plot_strain(
    compression_data, system, value_idx_at_targert_strain_new, time_at_target_strain_new, folder
):
    '''
    Plots the compression strain over time. 

    Args
    compression_data: pandas dataframe with columns Time_1e6, Time, Box-X, Compression strain
    system: int, the system number to analyze
    value_idx_at_targert_strain_new: float, the actual compression strain at the target strain
    time_at_target_strain_new: float, the time at which the the strain is the closest to the target strain
    folder: Path, the path to the specific compression folder


    Returns
    A plot the compression strain over time and save the figure. 
    The title of the plot has the target strain, the time at which the strain is closest to the target strain, and the actual strain at that time. 
    
    '''
    # combined plot time vs compression
    plt.plot(compression_data["Time_1e6"], compression_data["Compression strain"])
    plt.xlabel("time (ps)")
    plt.ylabel("compression strain")
    plt.title(
        f"system{system} target-strain{target_strain} real-strain {value_idx_at_targert_strain_new:.8g} at {time_at_target_strain_new}ps "
    )

    # Set major ticks every 0.1 on y-axis
    plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))

    # Draw horizontal red line at compression = 0.1
    plt.axhline(y=0.1, color="red", linestyle="--", linewidth=1)
    plt.axhline(y=0.2, color="orange", linestyle="--", linewidth=1)
    plt.axhline(y=0.6, color="green", linestyle="--", linewidth=1)

    plt.tight_layout()
    plt.savefig(
        folder
        / f"compression_{target_strain}strain_system{system}_{time_at_target_strain_new}ps.png"
    )
    plt.close()




####execute of script 
for system in systems:
    compression_data, compression_at_target_strain, time_at_target_strain, folder = ( 
        compression_plot_set_Lo(system)
    )

    compression_plot_strain(
        compression_data,
        system,
        compression_at_target_strain,
        time_at_target_strain,
        folder
    )
