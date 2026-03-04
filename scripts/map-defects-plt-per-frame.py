#plot defects
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
import matplotlib.colors as mcolors
import aggregate_defect_bar_plots as ag



configuration = "flat"
defect_cut_off = 10
restrictions = [200]

# lipid_compositions = {
#     "DOPC": [2, 0],
#     "DPPC": [2, 1],
#     "DOPC-DOPA": [1, 0],
#     "DPPC-DPPA": [1, 1],
#     "DOPC-DOPS": [0, 0],
#     "DPPC-DPPS": [0, 1],
# }

lipid_compositions = {
    "DOPC": [2, 0],
    "DPPC": [2, 1]
}

lipid_comp = list(lipid_compositions.keys())

#pulls out a frame that is representative of the average defect coverage 
defect_restricted_avg_system = ag.surface_coverage_defects(restrictions)
max_frames_avg_defect_coverage_per_system = ag.frame_extraction_avg_surface_defect_coverage(defect_restricted_avg_system)
print(max_frames_avg_defect_coverage_per_system)

#hardcoded paths
analysis_path = Path("/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-100ps-1200ext")
analysis_defect_path = analysis_path/f"defect-cut-off-{defect_cut_off}A"
compiled_data = Path(f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/defect-data-100ps-1200ext/{configuration}/defect-cut-off-{defect_cut_off}A/defect-maps")
compiled_data.mkdir(exist_ok=True)

# Loop through frames
# for number in range(5960,5961):
with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):

    for composition in lipid_comp:
        fig, ax_array = ph.fixed_size_subplots(1,1, subwidth=1.3, subheight=1,wmargin=0.5, hmargin=0.002, rmargin_scale = 0.1,tmargin_scale = 0.2)
        
        folder_system = analysis_path/f"defect-cut-off-{defect_cut_off}A/{composition}-{configuration}-production-stripped-ext"
        defect_counts = np.load(folder_system/"output-defects-upper.npy")
        frame_number = max_frames_avg_defect_coverage_per_system[composition]
        defect_counts_frame = defect_counts[frame_number]

        shape = np.load(folder_system/"output-shape.npy")
        X = (np.load(folder_system/"output-Xshape.npy"))/10
        Y = (np.load(folder_system/"output-Yshape.npy"))/10

        # center X so midpoint is at 0
        X_center = (X.max() + X.min()) / 2
        X = X - X_center

        reshaped_defects = defect_counts_frame.reshape(shape)

        masked_shallow_defects = (np.ma.masked_where(
            (reshaped_defects >= 1) | (reshaped_defects == 0),
            reshaped_defects,
        ))/10

        masked_deep_defects = (np.ma.masked_where(
            reshaped_defects != 0,
            reshaped_defects,
        ))/10

        
        # # get subplot position
        # i, j = position
        # ax = ax_array[i, j]  # pick the correct Axes

        # plot
        #cmap_single_orange = mcolors.ListedColormap(["orange"], name="constant_red")
        cmap_single_blue = mcolors.ListedColormap(["#1f91bc"], name="constant_blue")
        ax_array.pcolor(X, Y, masked_shallow_defects, cmap=cmap_single_blue, shading="auto")
        ax_array.pcolor(X, Y, masked_deep_defects, cmap=cmap_single_blue, shading="auto")

        # formatting
        ax_array.set_aspect("equal", adjustable="datalim")
        ax_array.xaxis.set_major_locator(MultipleLocator(5))
        ax_array.xaxis.set_minor_locator(MultipleLocator(2.5))
        ax_array.set_xlabel("X (nm)")
        ax_array.set_ylabel("Y (nm)")
        ax_array.set_xlim(X.min(), X.max())
        ax_array.set_ylim(Y.min(), Y.max())
        ax_array.set_aspect('equal', adjustable='box') 
        ax_array.set_title(composition, fontsize=16)  # just use ax, not ax[i,j]

    # save compiled frame
    
        plt.savefig(compiled_data/f"{composition}-{configuration}-{restrictions}-avg-defect-coverage-frame{frame_number}.pdf")
        plt.close(fig)
