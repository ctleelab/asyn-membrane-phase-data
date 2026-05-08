from pathlib import Path 
import MDAnalysis as mda 
import matplotlib.pyplot as plt
import numpy as np 
import ctleelab_plothelper.plothelpers as ph

configurations =  ["flat","strain2"]

start_box_dimX= {}
start_box_dimY= {}

lipid_colors = {
    "DOPC": "#1f77b4",
    "DPPC": "#ff7f0e",
    "DOPC-DOPA": "#2ca02c",
    "DPPC-DPPA": "#d62728",
    "DOPC-DOPS": "#9467bd",
    "DPPC-DPPS": "#8c564b"
}

for configuration in configurations:
    lipid_compositions = {
        "DOPC": 1 ,
        "DPPC": 2,
        "DOPC-DOPA": 6,
        "DPPC-DPPA": 5,
        "DOPC-DOPS": 3,
        "DPPC-DPPS": 4
    }

    for lipid, system in lipid_compositions.items(): 
        simulations_folder =  Path(f"/scratch/local/casakurai/asyn-phase-binding-data/simulations/curvature_selection/NVT")
        system_folder = simulations_folder/f"{configuration}"/f"system{system}-8x8x25-{lipid}-{configuration}-NVT"
        analysis_simulation_folder = simulations_folder/"analysis"
        analysis_simulation_folder.mkdir(exist_ok=True)

        if configuration == "flat":
            gro_file = system_folder/f"system{system}-large-flat-extractedgro.gro"


        if configuration == "strain2":
            gro = system_folder/"large-tessellation.gro"

        u = mda.Universe(gro_file)
        lx, ly, lz, alpha, beta, gamma = u.dimensions

        #covert to A
        lx = lx/10
        ly = ly/10
        lz = lz/10
        start_box_dimX[f"{lipid}-{configuration}"] = lx
        start_box_dimY[f"{lipid}-{configuration}"] = ly

#plot X0dim
systems = list(start_box_dimX.keys())
lx_values = list(start_box_dimX.values())

colors = [lipid_colors[system.split("-strain")[0].split("-flat")[0]] for system in systems]

colors = []
for system in systems:
    lipid = system.rsplit("-", 1)[0]   # removes "-flat" or "-strain2"
    colors.append(lipid_colors[lipid])

with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
    fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4.2, subheight=4.2)

    ax.bar(systems, lx_values,color = colors)

    ax.tick_params(axis="x", rotation=45)

    ax.set_xlabel("System")
    ax.set_ylabel("Starting Lx (nm)")

    fig.savefig(f"{analysis_simulation_folder}/Box-Lx-dim.png", dpi=300, bbox_inches="tight")

plt.close()

#plot Y-dim
systems = list(start_box_dimY.keys())
ly_values = list(start_box_dimY.values())

with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
    fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4.2, subheight=4.2)

    ax.bar(systems, ly_values, color = colors)

    ax.tick_params(axis="x", rotation=45)

    ax.set_xlabel("System")
    ax.set_ylabel("Starting Ly (nm)")

    fig.savefig(f"{analysis_simulation_folder}/Box-Ly-dim.png", dpi=300, bbox_inches="tight")

plt.close()

