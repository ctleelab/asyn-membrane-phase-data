

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
from scipy import ndimage as scipy_ndimage
from math import log
from collections import Counter
import pandas as pd
from scipy.stats import linregress


configurations = ["flat"]
#time = "100ps-1800ext"
defect_type = ["deep", "shallow", "all"]
lower_limit_size_defect = 15
upper_limit_size_defect = 200
restrictions = [200]
lipid_type = "DPPC-DPPS"
cut_off = "15A"
defect_cut_off = 6.3


#hardcoded paths
comparison_data = Path(f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/comparison-Packmem-methods")

lipid_compositions = {
    "DPPC-DPPS": [0, 0]
}




lipid_comp = list(lipid_compositions.keys())


def total_num_defect_restricted_X_cutoff(defect_data,shape,centroids,restriction,configuration,cutoff):
    
    defect_grid = defect_data.reshape(shape)

    #center of grids at 0 
    centered_centroids = centroids - np.mean(centroids)

    centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()

    centroid_mask_grid = centroid_X_mask.reshape(shape)

    restricted_region_grids_with_defects =  np.where((centroid_mask_grid  == 1) & (defect_grid  == 1),1,0)
    

    #here handle the defect based on size 
    all_defect_sizes_cutoff = []
    label_all, nfeat_all = scipy_ndimage.label(restricted_region_grids_with_defects)
    for i in range(1, nfeat_all + 1):
        defect_size = np.count_nonzero(label_all == i)
        if defect_size >= cutoff: 
            all_defect_sizes_cutoff.append(defect_size)


    total_grids = len(defect_data)

    defect_grids = np.sum(all_defect_sizes_cutoff)

    percent_defect_coverage = (defect_grids/total_grids)*100
    
    #return avg defect size 
    defect_avg_size = np.mean(all_defect_sizes_cutoff)
    #print("avg-defect-size:", defect_avg_size )




    #return total_num_defects_restricted_X, total_grids_restricted_X
    return defect_grids,  total_grids,  percent_defect_coverage, defect_avg_size 

# restrict the area to be analyzed
# then determine the size of the defects
# then only add the defects that are > cutoff
def surface_coverage_defects_cutoff(restrictions,cutoff): 
    defect_restricted_avg_system = {} 
    defect_size_restricted_avg_system = {}
    defect_restricted_per_frame_system = {}
    frames_avg_defect_coverage_per_system = {}
    for composition in lipid_comp:
        configurations_types = {}
        configurations_types_per_frame = {}
        configurations_types_frames_to_avg_defect = {}
        configurations_types_defect_size = {}

        for configuration in configurations: 
            if configuration ==  "flat": 
                time = "100ps-1800ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
            if configuration == "strain2":
                time = "100ps-1800ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"
            
            
            defect_counts = np.load(folder_system/"output-defects-upper.npy")
            quad_centroids_np = np.load(folder_system/"output-centroids.npy")
            shape = np.load(folder_system/"output-shape.npy")
            X = np.load(folder_system/"output-Xshape.npy")
            Y = np.load(folder_system/"output-Yshape.npy")

  
            n_frames = len(quad_centroids_np)

            percent_cut_defects = {}
            percent_cut_defects_per_frame  = {}
            avg_defect_size_over_frames = {}

            for restriction in restrictions: 

                total_defects_restricted_region_per_frame = np.empty(n_frames)
                total_grids_restricted_region_per_frame = np.empty(n_frames)
                percentage_defect_coverage_per_frame = np.empty(n_frames)
                avg_defect_size_per_frame = np.empty(n_frames)
                list_frames_close_avg_defect = []


                for frame in np.arange(0, n_frames):

                    defect_counts_frame = defect_counts[frame]
                    
                    centroids_frame =  quad_centroids_np[frame, :, 0:1]

                    # #all defect mask 
                    # masked_all_defects_num = np.where(
                    #     defect_counts_frame < 1,
                    #     1,
                    #     0
                    # )

                    masked_all_defects_num = np.where(defect_counts_frame < 1, 1, 0)


                    total_per_frame, total_grids, percent_defect_coverage, avg_defect_size = total_num_defect_restricted_X_cutoff(masked_all_defects_num,shape,centroids_frame,restriction,configuration,cutoff)

                    total_defects_restricted_region_per_frame[frame] = total_per_frame

                    total_grids_restricted_region_per_frame [frame] = total_grids

                    percentage_defect_coverage_per_frame[frame] = percent_defect_coverage

                    
                    # handles frames with no defects
                    if np.isnan(avg_defect_size):
                        avg_defect_size_per_frame[frame] = 0
                    else:
                        avg_defect_size_per_frame[frame] = avg_defect_size


                average_defect_coverage = np.mean(percentage_defect_coverage_per_frame)

    
                
                #print(f"{composition} avg defect coverage {average_defect_coverage}")
                std_defect_coverage = np.std(percentage_defect_coverage_per_frame)

                average_defect_size = np.mean(avg_defect_size_per_frame)
                average_defect_size_std = np.std(avg_defect_size_per_frame)

                shape = percentage_defect_coverage_per_frame.shape
                num_frames = shape[0]

                percent_cut_defects[restriction]  = [average_defect_coverage,std_defect_coverage,num_frames]
                avg_defect_size_over_frames[restriction] = [average_defect_size, average_defect_size_std ,num_frames]

                #storing per frame array, to look at convergence
                percent_cut_defects_per_frame[restriction] = percentage_defect_coverage_per_frame


            configurations_types[configuration] = percent_cut_defects
            configurations_types_defect_size[configuration] = avg_defect_size_over_frames

            configurations_types_per_frame[configuration] = percent_cut_defects_per_frame
            configurations_types_frames_to_avg_defect[configuration] = list_frames_close_avg_defect


        defect_restricted_avg_system[composition]= configurations_types
        defect_restricted_per_frame_system[composition] = configurations_types_per_frame
        frames_avg_defect_coverage_per_system[composition] = configurations_types_frames_to_avg_defect

        defect_size_restricted_avg_system[composition]= configurations_types_defect_size

    return defect_restricted_avg_system, defect_size_restricted_avg_system
defect_coverage_restricted_avg_system, defect_size_restricted_avg_system= surface_coverage_defects_cutoff(restrictions,15)




#####################################################################
###plot for the probability of finding a defect of a certain size### 
#####################################################################

#plotting new PackMem data 
major_probs = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]) 
restriction = 200
for composition, position in lipid_compositions.items():
    for configuration in configurations:
        if configuration == "flat":
            time = "100ps-1800ext"
            analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
            analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
            folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
        if configuration == "strain2":
            time = "100ps-1800ext"
            analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
            analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
            folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"

        # create a new figure for each system
        fig, ax = plt.subplots()
        defect_counts = np.load(folder_system/"output-defects-upper.npy")
        quad_centroids_np = np.load(folder_system/"output-centroids.npy")
        n_frames = len(quad_centroids_np)
        shape = np.load(folder_system/"output-shape.npy")

        defect_sizes = []
        bins = np.arange(0,4000,1)

        for frame in range(n_frames):

            defects_counts_frame = defect_counts[frame].reshape(shape)

            masked_all_defects_num = np.where(defects_counts_frame < 1, 1, 0)

            centroids_frame = quad_centroids_np[frame, :, 0:1]

            centered_centroids = centroids_frame - np.mean(centroids_frame)
            centroid_X_mask = np.where(
                (centered_centroids >= -restriction) & 
                (centered_centroids <= restriction), 1, 0
            ).squeeze()

            centroid_mask_grid = centroid_X_mask.reshape(shape)

            restricted_region_grids_with_defects = np.where(
                (centroid_mask_grid == 1) & 
                (masked_all_defects_num == 1), 1, 0
            )

            label_all, nfeat_all = scipy_ndimage.label(restricted_region_grids_with_defects)

            for i in range(1, nfeat_all + 1):
                defect_size = np.count_nonzero(label_all == i)
                if defect_size >= 15:
                    defect_sizes.append(defect_size)

        # histogram
        counts, _ = np.histogram(defect_sizes, bins=bins)
        probabilities = counts / counts.sum()
        bin_centers = 0.5 * (bins[1:] + bins[:-1])

        mask = probabilities > 0
        areas_nonzero = bin_centers[mask]
        probabilities_nonzero = probabilities[mask]
        log_probabilities = np.log(probabilities_nonzero)

        # apply PackMem-like fit limits
        masks_limits = (
            (areas_nonzero >= lower_limit_size_defect) &
            (probabilities_nonzero >= 1e-4)
        )


        areas_nonzero = areas_nonzero[masks_limits]
        probabilities_nonzero = probabilities_nonzero[masks_limits]
        log_probabilities = np.log(probabilities_nonzero)



        # scatter
        ax.scatter(areas_nonzero, log_probabilities, s=5, alpha=0.7, color = "black")

        # linear regression
        linear_regression = linregress(areas_nonzero, log_probabilities)

        slope = linear_regression.slope
        intercept = linear_regression.intercept
        r2_value = linear_regression.rvalue**2
        packdef_constant = abs(1 / slope)

        xfit = np.linspace(min(areas_nonzero), max(areas_nonzero), 100)
        yfit = slope * xfit + intercept

        ax.plot(xfit, yfit, color="black",
                label=f"New-code\ny={slope:.2f}x+{intercept:.2f}\n$R^2$={r2_value:.2f}\n defect-constant={packdef_constant:.2f}")

        # axis formatting
        ax.set_yticks(np.log(major_probs))
        ax.set_yticklabels([
            r"$10^{-1}$", r"$10^{-2}$", r"$10^{-3}$",
            r"$10^{-4}$", r"$10^{-5}$", r"$10^{-6}$"
        ])

        ax.set_ylabel("ln(Probability)")
        ax.set_xlabel("Defect Area (Å²)")
        ax.set_ylim(np.log(1e-6), 0)


        # plt.savefig(
        #     f"{compiled_data}/{composition}-{configuration}-defect-size-probability.png"
        # )

    #plot old packMem data
    Packmem_old_path = Path("/scratch/local/casakurai/asyn-phase-binding-data/analysis/Packmem-comparison-old")
    Packmem_old_data = f"{Packmem_old_path}/{composition}-5fr-Total-Up-all.txt"

    #read in data 
    import numpy as np

    packdef_data_size = np.loadtxt(Packmem_old_data, comments="#")[:,1]


    # histogram
    bins = np.arange(0,4000,1)
    counts, _ = np.histogram(packdef_data_size, bins=bins)
    probabilities = counts / counts.sum()
    bin_centers = 0.5 * (bins[1:] + bins[:-1])

    mask = probabilities > 0
    areas_nonzero = bin_centers[mask]
    probabilities_nonzero = probabilities[mask]
    log_probabilities = np.log(probabilities_nonzero)

    # apply PackMem-like fit limits
    masks_limits = (
            (areas_nonzero >= lower_limit_size_defect) &
            (probabilities_nonzero >= 1e-4)
        )


    areas_nonzero = areas_nonzero[masks_limits]
    probabilities_nonzero = probabilities_nonzero[masks_limits]
    log_probabilities = np.log(probabilities_nonzero)


     

    # scatter
    ax.scatter(areas_nonzero, log_probabilities, s=5, alpha=0.7,color = "red")

    # linear regression
    linear_regression = linregress(areas_nonzero, log_probabilities)

    slope = linear_regression.slope
    intercept = linear_regression.intercept
    r2_value = linear_regression.rvalue**2
    packdef_constant = abs(1 / slope)

    xfit = np.linspace(min(areas_nonzero), max(areas_nonzero), 100)
    yfit = slope * xfit + intercept

    ax.plot(xfit, yfit, color="red",
            label=f"Old-code\ny={slope:.2f}x+{intercept:.2f}\n$R^2$={r2_value:.2f}\n defect-constant={packdef_constant:.2f}")

    ax.legend(fontsize="small")
    ax.set_title(f"{composition} — {configuration}")

    plt.tight_layout()


    plt.savefig(
    f"{comparison_data}/{composition}-5fr-comparison-defect-size-probability.png"
    )
    
    plt.close()
