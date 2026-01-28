#plot defects
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
import matplotlib.colors as mcolors



configuration = "flat"
defect_cut_off = 15

lipid_compositions = {
    "DOPC": [2, 0],
    "DPPC": [2, 1],
    "DOPC-DOPA": [1, 0],
    "DPPC-DPPA": [1, 1],
    "DOPC-DOPS": [0, 0],
    "DPPC-DPPS": [0, 1],
}


analysis_path = Path("/home/casakurai/scratch/asyn-phase-binding-data/analysis/defect-data")
analysis_defect_path = analysis_path/f"defect-cut-off-{defect_cut_off}A"
compiled_data = Path(f"/home/casakurai/scratch/asyn-phase-binding-data/Figures/defect-data/{configuration}/defect-cut-off-{defect_cut_off}A")
compiled_data.mkdir(exist_ok=True)

# Loop through frames
for number in range(0,31):
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax_array = ph.fixed_size_subplots(3, 2, subwidth=3, subheight=1.5,wmargin=0.5, hmargin=0.002, rmargin_scale = 0.1,tmargin_scale = 0.2)

        for composition, position in lipid_compositions.items():
            folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
            defect_counts = np.load(folder_system/"output-defects-upper.npy")
            defect_counts_frame = defect_counts[number]

            shape = np.load(folder_system/"output-shape.npy")
            X = np.load(folder_system/"output-Xshape.npy")
            Y = np.load(folder_system/"output-Yshape.npy")

            reshaped_defects = defect_counts_frame.reshape(shape)

            masked_shallow_defects = np.ma.masked_where(
                (reshaped_defects >= 1) | (reshaped_defects == 0),
                reshaped_defects,
            )

            masked_deep_defects = np.ma.masked_where(
                reshaped_defects != 0,
                reshaped_defects,
            )

            # get subplot position
            i, j = position
            ax = ax_array[i, j]  # pick the correct Axes

            # plot
            cmap_single_orange = mcolors.ListedColormap(["orange"], name="constant_orange")
            ax.pcolor(X, Y, masked_shallow_defects, cmap=cmap_single_orange, shading="auto")

            cmap_single_red = mcolors.ListedColormap(["red"], name="constant_red")
            ax.pcolor(X, Y, masked_deep_defects, cmap=cmap_single_red, shading="auto")

            # formatting
            ax.set_aspect("equal", adjustable="datalim")
            ax.xaxis.set_major_locator(MultipleLocator(25))
            ax.set_xlabel("X (A)")
            ax.set_ylabel("Y (A)")
            ax.set_xlim(X.min(), X.max())
            ax.set_ylim(Y.min(), Y.max())
            ax.set_aspect('equal', adjustable='box') 
            ax.set_title(composition, fontsize=16)  # just use ax, not ax[i,j]

        # save compiled frame
        
        plt.savefig(compiled_data/f"frame{number}")
        print(f"finished frame {number}")
    plt.close(fig)
