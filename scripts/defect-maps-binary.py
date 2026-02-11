
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator



configuration = "flat" 
defect_cut_off = 10
time = "100ps"

defect_type = ["deep", "shallow", "all"]

analysis_path = f"/home/casakurai/scratch/asyn-phase-binding-data/analysis/defect-data-{time}"
# analysis_path = Path(analysis_path)
analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")



lipid_compositions = {
    "DOPC": [0, 0],
    "DPPC": [0, 1],
    "DOPC-DOPA": [1, 0],
    "DPPC-DPPA": [1, 1],
    "DOPC-DOPS": [2, 0],
    "DPPC-DPPS": [2, 1]
}



for composition in lipid_compositions:
    folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped"
    defect_counts = np.load(folder_system/"output-defects-upper.npy")
    quad_centroids_np = np.load(folder_system/"output-centroids.npy")
    shape = np.load(folder_system/"output-shape.npy")
    X = np.load(folder_system/"output-Xshape.npy")
    Y = np.load(folder_system/"output-Yshape.npy")

    frame =  0


    fig, ax = ph.fixed_size_subplots(
        1, 1, subwidth=1.3, subheight=.5, rmargin_scale=1.5
    )

    defects_frame = defect_counts[frame]

    reshaped_defects = defects_frame.reshape(shape)

    masked_shallow_defects = np.ma.masked_where(
        (reshaped_defects >= 1) | (reshaped_defects == 0),
        reshaped_defects,
    )
    cmap_single_orange = mcolors.ListedColormap(
        ["orange"], name="constant_orange"
    )

    X_centered = X - X.mean()
    X_centered_nm = X_centered/10
    Y_nm = Y/10
    
    ax.pcolor(
         X_centered_nm, Y_nm, masked_shallow_defects, cmap=cmap_single_orange, shading="auto"
    )

    masked_deep_defects = np.ma.masked_where(
        reshaped_defects != 0,
        reshaped_defects,
    )
    cmap_single_red = mcolors.ListedColormap(["red"], name="constant_red")
    
    ax.pcolor(X_centered_nm, Y_nm, masked_deep_defects, cmap=cmap_single_red, shading="auto")
    
    ax.set_aspect("equal", adjustable="datalim")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax.set_xlabel("X (nm)")
    ax.set_ylabel("Y (nm)")
    fig.savefig(f"{frame}frame-{composition}-{configuration}-upper-defect.pdf")

    plt.close(fig)