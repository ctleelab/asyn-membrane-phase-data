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


systems = [1]
initial_dim = "8x8x25"
final_dim = "large"
sim_time = 800
pressure = "3"
step = 7.2
compression = "35bar"
frame = 85620
ensemble = "NPT"
lipid = "DPPC"

for sys in systems:
    #file set-up when files are in system folder
    system_folder = f"system{sys}-{initial_dim}"
    system_path = base_path/system_folder
    # compression_folder = f"xzPcoupled-{pressure}bar-compression-{sim_time}ns"
    #NVT compression_folder
    compression_folder = f"xzPcoupled-{pressure}bar-compression-{sim_time}ns"
    compression_path = system_path / compression_folder
    edr = compression_path /"large-compression.edr"
    analysis_folder= analysis_path/f"testing_compression"
    analysis_folder.mkdir(exist_ok=True)

    # #file set-up when files are in analysis folder
    # analysis = util.analysis_path
    # analysis_curv_folder = f"curvature_selection"
    # system_analysis_folder = analysis/analysis_curv_folder/f"{compression}_{sim_time}ns"/f"system{sys}-{initial_dim}-{lipid}"
    # system_analysis_folder.mkdir(exist_ok = True)
    # system_analysis_time_folder = f"system{sys}-{initial_dim}-{lipid}-{frame}ps-{ensemble}"
    # analysis_curv_dir = analysis/analysis_curv_folder/f"{compression}_{sim_time}ns"/f"{system_analysis_folder}"/f"{system_analysis_time_folder}"
    # analysis_curv_dir.mkdir(exist_ok = True)
    # edr = analysis_curv_dir /"equil"/f"equilibration{step}.edr"


    aux = mda.auxiliary.EDR.EDRReader(edr)
    terms = aux.get_data(["Box-X", "Box-Z", "Pressure", "Pres-XX", "Pres-YY", "Pres-ZZ"])
    #terms = aux.get_data(["Pressure", "Pres-XX", "Pres-YY", "Pres-ZZ"])


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
    #df = pd.DataFrame({"Time": terms["Time"], "Pressure": terms["Pressure"], "Pres-XX": terms["Pres-XX"], "Pres-YY": terms["Pres-YY"], "Pres-ZZ": terms["Pres-ZZ"]})
    #df.to_csv(analysis_folder/f"compression_strain_system{sys}_{final_dim}{initial_dim}_{pressure}bar_{sim_time}ns.csv", index = False)

    target_strain = .6
    #find the index where the compression strain is closest to .2
    strain_idx = np.argmin(np.abs(df["Compression strain"] - target_strain))
    time_at_target_strain = df["Time"].iloc[strain_idx]
    print(f"The strain is {target_strain} at {time_at_target_strain}ps")


    #  #combined plot of box-X vs box-Z dimension 
    # plt.plot(df["Time"], df["Box-X"])
    # plt.plot(df["Time"],df["Box-Z"])
    # plt.xlabel('Time')
    # plt.ylabel('dimension (Nm)')
    # plt.legend(['x dimension', 'z dimension'])
    # plt.title(f"box x and z dimension for system{sys} at {pressure}bar")
    # plt.savefig(analysis_folder / f"large-xzdimension{sys}-{pressure}bar-xz-Pcoupled-{sim_time}ns.png")
    # plt.close()


    #combined plot time vs compression
    plt.plot(df["Time_1e6"], df["Compression strain"])
    plt.xlabel('time (ps)')
    plt.ylabel('compression strain')
    plt.title(f"compression strain for system{sys} {final_dim}{initial_dim} at {pressure}bar for 800ns")

    # Set major ticks every 0.1 on y-axis
    plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))

    # Draw horizontal red line at compression = 0.1
    plt.axhline(y=0.1, color='red', linestyle='--', linewidth=1)
    plt.axhline(y=0.2, color='orange', linestyle='--', linewidth=1)
    plt.axhline(y=0.6, color='green', linestyle='--', linewidth=1)


    plt.tight_layout()
    plt.savefig(analysis_folder / f"compression_strain_system{sys}_{final_dim}{initial_dim}_{pressure}bar_800ns.png")
    #plt.savefig(analysis_folder / f"3bar_160104ps_imposed_compression_strain_system{sys}_{final_dim}{initial_dim}_{pressure}bar_1000ns.png")
    plt.close()

    # # #combined plot time vs pressure 
    # plt.plot(df["Time"], df["Pressure"])
    # plt.xlabel('time(ps)')
    # plt.ylabel('Pressure')
    # #plt.title(f"Pressure for system{sys} {final_dim}{initial_dim} at NVT for 500ns")
    # plt.title(f"Pressure for system{sys} {final_dim}{initial_dim} at NVT compression {compression}")
    # #plt.savefig(analysis_folder / f"Pressure_system{sys}_{final_dim}{initial_dim}_bar_500ns.png")
    # figures = analysis/analysis_curv_folder/f"{compression}_{sim_time}ns"/"Pressure_plot"
    # plt.savefig(figures/ f"3bar_{frame}ps_imposed_compression_strain_system{sys}_{final_dim}{initial_dim}_{pressure}bar_{ensemble}-equilibration{step}.png")
    # plt.close()

    # #length-X vs compression 
    # plt.plot(df["Box-X-nm"], df["Compression strain"])
    # plt.xlabel('Box-X (nm)')
    # plt.ylabel('compression strain')
    # plt.title(f"compression strain for system{sys} {final_dim}{initial_dim} at {pressure}bar for {sim_time}ns")
    # plt.savefig(analysis_folder / f"compression_strain_vs_Xdim_system{sys}_{final_dim}{initial_dim}_{pressure}bar_{sim_time}ns.png")
    # plt.close()

    

    # #time vs X and Z length 
    # #use to determine if the box size changes, when I go from DOPC to DPPC
    # plt.plot(df["Time"], df["Box-Z-nm"])
    # plt.plot(df["Time"], df["Box-X-nm"])
    # plt.xlabel('time(ps)')
    # plt.ylabel('Box dimension (nm)')
    # plt.legend(["Z-dim", "X-dim"])
    # plt.title(f"compression strain for system{sys} {final_dim}{initial_dim} at {pressure}bar for {sim_time}ns")
    # plt.savefig(analysis_folder / f"time_vs_ZandXdim_system{sys}_{final_dim}{initial_dim}_{pressure}bar_{sim_time}ns.png")
    # plt.close()