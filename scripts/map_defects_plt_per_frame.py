# plot defects
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
import matplotlib.colors as mcolors
import defect_coverage as dc
import pandas as pd

configurations = ["flat","strain2"]
packMem_cutoff = 6.3
restriction = 200
defect_size_cut_off = 15
PackMem_frame_interval = 5 

lipid_compositions = {
    "DOPC": [2, 0],
    "DPPC": [2, 1],
    "DOPC-DOPA": [1, 0],
    "DPPC-DPPA": [1, 1],
    "DOPC-DOPS": [0, 0],
    "DPPC-DPPS": [0, 1],
}


lipid_comp = list(lipid_compositions.keys())

figures_data = {}

# create defect map plots
with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):

    for configuration in configurations:
        # hardcoded paths
        analysis_path = Path(
            "/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-100ps-1800ext"
        )


        analysis_defect_path = analysis_path / f"defect-cut-off-{packMem_cutoff}A"
        configuration_data = Path(
            f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/defect-data-100ps-1800ext/{configuration}/defect-cut-off-{packMem_cutoff}A/defect-maps"
        )
        configuration_data.mkdir(exist_ok=True)

        compiled_data = Path(
            f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/defect-data-100ps-1800ext/compiled/defect-cut-off-{packMem_cutoff}A/defect-maps"
        )
        compiled_data.mkdir(exist_ok=True)

        for composition in lipid_comp:
            
            #load the defect data per system
            if configuration == "flat":
                folder_system = (
                    analysis_path
                    / f"defect-cut-off-{packMem_cutoff}A/{composition}-{configuration}-production-stripped-ext1800-cutoff"
                )
            if configuration == "strain2":
                folder_system = (
                    analysis_path
                    / f"defect-cut-off-{packMem_cutoff}A/{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"
                )

            # computes the average defect coverage over the entire simulation
            (
                avg_percent_coverage_defects,
                avg_defect_size_over_frames,
                percent_coverage_defects_per_frame,
                n_frames,
            ) = dc.surface_coverage_defects_cutoff(
                restriction,
                defect_size_cut_off,
                packMem_cutoff,
                composition,
                configuration,
            )

            # determines the last frame of the trajectory that has the average defect coverage
            max_frames_avg_defect_coverage_per_system = (
                dc.frame_extraction_avg_surface_defect_coverage(
                    avg_percent_coverage_defects,
                    n_frames,
                    percent_coverage_defects_per_frame,
                )
            )

            fig, ax_array = ph.fixed_size_subplots(
                1,
                1,
                subwidth=1.3,
                subheight=1,
                wmargin=0.5,
                hmargin=0.002,
                rmargin_scale=0.1,
                tmargin_scale=0.2,
            )

            # load data
            defect_counts = np.load(folder_system / "output-defects-upper.npy")
    
            defect_counts_frame = defect_counts[max_frames_avg_defect_coverage_per_system]

            shape = np.load(folder_system / "output-shape.npy")
            X = (np.load(folder_system / "output-Xshape.npy")) / 10
            Y = (np.load(folder_system / "output-Yshape.npy")) / 10

            # center X so midpoint is at 0
            X_center = (X.max() + X.min()) / 2
            X = X - X_center

            reshaped_defects = defect_counts_frame.reshape(shape)

            # mask for deep and shallow defects
            masked_shallow_defects = (
                np.ma.masked_where(
                    (reshaped_defects >= 1) | (reshaped_defects == 0),
                    reshaped_defects,
                )
            ) / 10

            masked_deep_defects = (
                np.ma.masked_where(
                    reshaped_defects != 0,
                    reshaped_defects,
                )
            ) / 10

            # plot
            # cmap_single_orange = mcolors.ListedColormap(["orange"], name="constant_red")
            cmap_single_blue = mcolors.ListedColormap(["#1f91bc"], name="constant_blue")
            ax_array.pcolor(
                X, Y, masked_shallow_defects, cmap=cmap_single_blue, shading="auto"
            )
            ax_array.pcolor(
                X, Y, masked_deep_defects, cmap=cmap_single_blue, shading="auto"
            )

            # formatting
            ax_array.set_aspect("equal", adjustable="datalim")
            ax_array.xaxis.set_major_locator(MultipleLocator(5))
            ax_array.xaxis.set_minor_locator(MultipleLocator(2.5))
            ax_array.set_xlabel("X (nm)")
            ax_array.set_ylabel("Y (nm)")
            ax_array.set_xlim(X.min(), X.max())
            ax_array.set_ylim(Y.min(), Y.max())
            ax_array.set_aspect("equal", adjustable="box")
            ax_array.set_title(composition, fontsize=16)  # just use ax, not ax[i,j]

            # save compiled frame
            plt.savefig(
                configuration_data
                / f"{composition}-{configuration}-{restriction}-avg-defect-coverage-frame{max_frames_avg_defect_coverage_per_system}.pdf"
            )
            plt.savefig(
                configuration_data
                / f"{composition}-{configuration}-{restriction}-avg-defect-coverage-frame{max_frames_avg_defect_coverage_per_system}.png"
            )
            plt.close(fig)



            print("max-frame:", max_frames_avg_defect_coverage_per_system)
            print("average defect coverage:", avg_percent_coverage_defects[0])
            MD_frame = max_frames_avg_defect_coverage_per_system*PackMem_frame_interval
            print("corresponding MD frame:",  MD_frame)

            lipid_configuration = f"{composition}-{configuration}"
            figures_data[lipid_configuration] = [max_frames_avg_defect_coverage_per_system, avg_percent_coverage_defects[0], MD_frame]


#export out figures_data as .csv 
rows = []
for key in figures_data:
    frame, average_percent_coverage_defect, MD_image = figures_data[key]
    rows.append({
        "system": key,
        "representative_defect_map_frame": frame,
        "avg_defect_coverage": average_percent_coverage_defect, 
        "MD_frame": MD_image
    })
df = pd.DataFrame(rows)
df.to_csv(f"{compiled_data}/reshaped_frames_visualization.csv", index = False)
