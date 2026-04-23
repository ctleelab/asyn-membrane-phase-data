

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


configurations = ["flat","strain.2"]
#configurations = ["strain.2"]
defect_cut_off = 6.3
time = "100ps-1800ext"
cutoff = 15 
defect_type = ["deep", "shallow", "all"]
restrictions = [25]
lipid_type = "all-systems"
cut_off = "15A" 

#hardcoded paths
compiled_data = Path(f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/defect-data-100ps-1800ext-equil-ext/compiled/defect-cut-off-{defect_cut_off}A")
#compiled_data_reshaped = Path(f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/defect-data-100ps-1800ext-equil-ext/compiled/defect-cut-off-{defect_cut_off}A/reshaped")
compiled_data_reshaped = Path(f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/trials")
compiled_data.mkdir(exist_ok=True)
compiled_data_reshaped.mkdir(exist_ok=True)

lipid_compositions = {
    "DOPC": [0, 0],
    "DPPC": [0, 1],
    "DOPC-DOPA": [1, 0],
    "DPPC-DPPA": [1, 1],
    "DOPC-DOPS": [2, 0],
    "DPPC-DPPS": [2, 1]
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
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}-equil-ext"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
            if configuration == "strain.2":
                time = "100ps-1800ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}-equil-ext"
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
                        avg_defect_size_per_frame[frame] = avg_defect_size/100 #convert to nm^2


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

    return defect_restricted_avg_system, defect_size_restricted_avg_system,defect_restricted_per_frame_system

defect_coverage_restricted_avg_system, defect_size_restricted_avg_system, defect_restricted_per_frame_system = surface_coverage_defects_cutoff(restrictions,cutoff)



def frame_extraction_avg_surface_defect_coverage(defect_restricted_avg_system,configuration,composition,defect_size_cutoff):
    '''
    Compiles a list of frames that is +/- .01 away from the avg defect coverage over the entire membrane

    Args:

    Returns:
    max_frames_avg_defect_coverage_per_system (int): A single frame that represents the avg coverage of defects across the entire membrane for each lipid composition.
   
    '''
    frames_avg_defect_coverage_per_system = {} 
    max_frames_avg_defect_coverage_per_system = {}
    for composition in lipid_comp:
        configurations_types_frames_to_avg_defect = {}
        list_frames_close_avg_defect = []

        for configuration in configurations:
            if configuration ==  "flat": 
                time = "100ps-1800ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}-equil-ext"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
            if configuration == "strain.2":
                time = "100ps-1800ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}-equil-ext"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"

            defect_counts = np.load(folder_system/"output-defects-upper.npy")
            quad_centroids_np = np.load(folder_system/"output-centroids.npy")
            shape = np.load(folder_system/"output-shape.npy")
            n_frames = len(quad_centroids_np)

            for restriction in restrictions: 
                for frame in np.arange(0, n_frames):
                    defect_counts_frame = defect_counts[frame]
                    centroids_frame =  quad_centroids_np[frame, :, 0:1]

                    #all defect 
                    masked_all_defects_num = np.where(
                        defect_counts_frame < 1,
                        1,
                        0
                    )

                #determines the defect percent coverage for frame
                    _, _,  percent_defect_coverage, _ = total_num_defect_restricted_X_cutoff(masked_all_defects_num,shape,centroids_frame,restriction,configuration,defect_size_cutoff)


                    if composition == "DOPC":
                        avg_defect_DOPC = defect_restricted_avg_system["DOPC"][configuration][restriction][0]
                        upper_bound_DOPC = avg_defect_DOPC + .01
                        lower_bound_DOPC = avg_defect_DOPC - .01
                        if lower_bound_DOPC < percent_defect_coverage < upper_bound_DOPC: #storing frames close to the average defect coverage for entire membrane 
                            list_frames_close_avg_defect.append(frame)
                    if composition == "DPPC":
                        avg_defect_DPPC = defect_restricted_avg_system["DPPC"][configuration][restriction][0]
                        upper_bound_DPPC = avg_defect_DPPC + .01
                        lower_bound_DPPC = avg_defect_DPPC - .01
                        if lower_bound_DPPC < percent_defect_coverage < upper_bound_DPPC: #storing frames close to the average defect coverage for entire membrane 6
                            list_frames_close_avg_defect.append(frame)
                    if composition == "DOPC-DOPA":
                        avg_defect_DOPC_DOPA = defect_restricted_avg_system["DOPC-DOPA"][configuration][restriction][0]
                        upper_bound_DOPC_DOPA = avg_defect_DOPC_DOPA + .01
                        lower_bound_DOPC_DOPA = avg_defect_DOPC_DOPA - .01
                        if lower_bound_DOPC_DOPA < percent_defect_coverage < upper_bound_DOPC_DOPA: #storing frames close to the average defect coverage for entire membrane 
                            list_frames_close_avg_defect.append(frame)
                    if composition == "DPPC-DPPA":
                        avg_defect_DPPC_DPPA = defect_restricted_avg_system["DPPC-DPPA"][configuration][restriction][0]
                        upper_bound_DPPC_DPPA  = avg_defect_DPPC_DPPA  + .01
                        lower_bound_DPPC_DPPA  = avg_defect_DPPC_DPPA  - .01
                        if lower_bound_DPPC_DPPA < percent_defect_coverage < upper_bound_DPPC_DPPA: #storing frames close to the average defect coverage for entire membrane 
                            list_frames_close_avg_defect.append(frame)
                    if composition == "DOPC-DOPS":
                        avg_defect_DOPC_DOPS = defect_restricted_avg_system["DOPC-DOPS"][configuration][restriction][0]
                        upper_bound_DOPC_DOPS = avg_defect_DOPC_DOPS + .01
                        lower_bound_DOPC_DOPS = avg_defect_DOPC_DOPS - .01
                        if lower_bound_DOPC_DOPS < percent_defect_coverage < upper_bound_DOPC_DOPS: #storing frames close to the average defect coverage for entire membrane 
                            list_frames_close_avg_defect.append(frame)
                    if composition == "DPPC-DPPS":
                        avg_defect_DPPC_DPPS = defect_restricted_avg_system["DPPC-DPPS"][configuration][restriction][0]
                        upper_bound_DPPC_DPPS = avg_defect_DPPC_DPPS  + .01
                        lower_bound_DPPC_DPPS  = avg_defect_DPPC_DPPS  - .01
                        if lower_bound_DPPC_DPPS < percent_defect_coverage < upper_bound_DPPC_DPPS: #storing frames close to the average defect coverage for entire membrane 
                            list_frames_close_avg_defect.append(frame)


        configurations_types_frames_to_avg_defect[configuration] = list_frames_close_avg_defect
        frames_avg_defect_coverage_per_system[composition] = configurations_types_frames_to_avg_defect
        max_frames_avg_defect_coverage_per_system[composition] = max(frames_avg_defect_coverage_per_system[composition][configuration]) 
    return max_frames_avg_defect_coverage_per_system



#######################################
#bar graphs for avg defect coverage ###
#######################################

def bar_graph_system_peak_grouped_defect_coverage(data,restriction,configurations):
    systems = list(data.keys())
    n_configs = len(configurations)

    #set labels for system
    x = np.arange(len(systems))  

    width = 0.35 
    system_spacing = 4   # separates systems
    config_spacing = 1  # keeps configs tight
    
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.2, subheight=1.5)
        

        #plot each configuration
        for i, config in enumerate(configurations):
            avg_values = [data[sys][i][0] for sys in systems]
            std_values = [data[sys][i][1] for sys in systems]
            group_spacing = 1.2 
            offsets = (i - (n_configs - 1) / 2) * width * config_spacing
            # print(systems)
            # print(f"{config}-avg-defect-coverage:", avg_values)

            ax.bar(x + offsets,  avg_values, width ,yerr = std_values, label=config, capsize=5)
        ax.legend()
        ax.set_ylabel('Surface Covered by Defects (%)')
        ax.set_ylim(0,20)
        ax.set_xticks(x)
        #ax.set_ylim(0,2)

        ax.set_xticklabels(lipid_comp,  ha='right')
        ax.tick_params(axis='x', rotation=45)

    plt.savefig(f"{compiled_data_reshaped}/reshaped-{config}-peak{restriction}-{lipid_type}-cutoff{cut_off}-percentcoverage.png")
    plt.savefig(f"{compiled_data_reshaped}/reshaped-{config}-peak{restriction}-{lipid_type}-cutoff{cut_off}-percentcoverage.pdf")

###################################
#bar graphs for avg defect size ###
###################################

def bar_graph_system_peak_grouped_defect_size(data,restriction,configurations):
    systems = list(data.keys())
    n_configs = len(configurations)

    #set labels for system
    x = np.arange(len(systems))  

    width = 0.35 
    config_spacing = 1  # keeps configs tight
    
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.2, subheight=1.5)
        

        #plot each configuration
        for i, config in enumerate(configurations):
            avg_values = [data[sys][i][0]  for sys in systems]
            std_values = [data[sys][i][1]  for sys in systems] 
            offsets = (i - (n_configs - 1) / 2) * width * config_spacing
            # print(systems)
            # print(f"{config}-avg-defect-coverage:", avg_values)

            ax.bar(x + offsets,  avg_values, width ,yerr = std_values, label=config, capsize=5)
        ax.legend()
        ax.set_ylabel('Average Defect Area (nm²)')
        ax.set_ylim(0,3)
        ax.set_xticks(x)
        #ax.set_ylim(0,2.3)

        ax.set_xticklabels(lipid_comp,  ha='right')
        ax.tick_params(axis='x', rotation=45)

        plt.savefig(f"{compiled_data_reshaped }/reshaped-{config}-peak{restriction}-{lipid_type}-cutoff{cut_off}-defectsize.png")
        plt.savefig(f"{compiled_data_reshaped }/reshaped-{config}-peak{restriction}-{lipid_type}-cutoff{cut_off}-defectsize.pdf")

# #######################################################
# ##execution code for bar graphs for avg defect size ###
# #######################################################
restriction = 25
data_grouped = {}
for system in lipid_comp:
    data_grouped[system] = []
    for config in configurations:
        avg, std, num_frames = defect_size_restricted_avg_system[system][config][restriction]
        data_grouped[system].append([avg, std])
bar_graph_system_peak_grouped_defect_size(data_grouped, restriction, configurations)

rows = []
for key in defect_size_restricted_avg_system:
    for config in configurations:
        avg, std, num_frames = defect_size_restricted_avg_system[key][config][restriction]
        variable_name = f"{key}_{config}_{restriction}"
        rows.append({
            "system": variable_name,
            "avg_defect_size(nm²)": avg,
            "std": std,
            "num_frames": num_frames
        })
df = pd.DataFrame(rows)
df.to_csv(f"{compiled_data_reshaped}/reshaped_{lipid_type}_defect_data_avg_size{restriction}-cutoff{cut_off}.csv", index = False)

# #######################################################
# ##execution code for bar graphs for defect coverage ###
# #######################################################
restriction = 25
data_grouped = {}
for system in lipid_comp:
    data_grouped[system] = []
    for config in configurations:
        avg, std, num_frames = defect_coverage_restricted_avg_system[system][config][restriction]
        data_grouped[system].append([avg, std])
bar_graph_system_peak_grouped_defect_coverage(data_grouped, restriction, configurations)
rows = []
for key in defect_coverage_restricted_avg_system:
    for config in configurations:
        avg, std, num_frames = defect_coverage_restricted_avg_system[key][config][restriction]
        variable_name = f"{key}_{config}_{restriction}"
        rows.append({
            "system": variable_name,
            "avg_defect_coverage(%)": avg,
            "std": std,
            "num_frames": num_frames
        })
df = pd.DataFrame(rows)
df.to_csv(f"{compiled_data_reshaped }/reshaped_{lipid_type}_defect_data_avg_coverage{restriction}_cutoff{cut_off}.csv", index = False)


# # # ################################################
# # # #scatter plot for defect coverage over frames###
# # # ################################################


def defect_coverage_simulation(defect_restricted_per_frame_system,restrictions, configurations):
    '''
    Plots the percent defect coverage for each frame for a given configuration and restriction, to look at convergence of defect coverage over time.
    Can handle two compositions at a time. 
    '''
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        for restriction in restrictions:
            for configuration in configurations:
                fig, ax = ph.fixed_size_subplots(1, 1, subwidth=5, subheight=3)

                #colors_presets = [ "#f8c362ff",  "#5db7fc" ] 


                for i, composition in enumerate(lipid_compositions):

                    per_frame_array = \
                        defect_restricted_per_frame_system[composition][configuration][restriction]

                    frames = np.arange(len(per_frame_array))

                    ax.scatter(
                        frames,
                        per_frame_array,
                        s=3,
                        label=composition
                    )
                    ax.set_ylim(0,22)
                    #ax.set_xlim(250,)
                    ax.set_xlabel("Frame")
                    ax.set_ylabel("Percent defect coverage")
                    ax.set_title(f"Defect coverage per frame ({configuration}, restriction={restriction} Å)")

                    ax.legend(frameon=False)

                    plt.tight_layout()

                    plt.savefig(
                        f"{compiled_data_reshaped }/frame-defect-coverage-{configuration}-restriction{restriction}-1800ext-{lipid_type}.png",
                        dpi=300
                    )
defect_coverage_simulation(defect_restricted_per_frame_system,restrictions,configurations)



#################################################################
###plot the probability of finding a defect of a certain size### 
##################################################################
major_probs = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]) 
restriction = 25
system_defect_sizes_log_transformed_prob = {}
for composition, position in lipid_compositions.items():
    strain_defect_sizes_log_transformed_prob = {}
    for configuration in configurations:
        if configuration ==  "flat": 
            time = "100ps-1800ext"
            analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}-equil-ext"
            analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
            folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
        if configuration == "strain.2":
            time = "100ps-1800ext"
            analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}-equil-ext"
            analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
            folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"
        
        defect_counts = np.load(folder_system/"output-defects-upper.npy")
        quad_centroids_np = np.load(folder_system/"output-centroids.npy")
        n_frames = len(quad_centroids_np)
        shape = np.load(folder_system/"output-shape.npy")

        defect_sizes = []
        bins = np.arange(0,2500,1)  #each bin is 1A^2

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
                (areas_nonzero >= cutoff) &
                (probabilities_nonzero >= 1e-4)
            )


            areas_nonzero = areas_nonzero[masks_limits]
            probabilities_nonzero = probabilities_nonzero[masks_limits]
            log_probabilities = np.log(probabilities_nonzero)


            #saves the log probability data per configuration 
            strain_defect_sizes_log_transformed_prob[configuration]= [areas_nonzero,log_probabilities]
        
        system_defect_sizes_log_transformed_prob[composition] = strain_defect_sizes_log_transformed_prob

