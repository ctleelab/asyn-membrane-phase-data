

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


configurations = ["flat","strain2"]
defect_cut_off = 6.3
time = "100ps-1800ext"
defect_type = ["deep", "shallow", "all"]
lower_limit_size_defect = 15
restrictions = [40]
lipid_type = "all-systems"
cut_off = "15A"


#hardcoded paths
compiled_data = Path(f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/defect-data-100ps-1800ext/compiled/defect-cut-off-{defect_cut_off}A")
compiled_data_reshaped = Path(f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/defect-data-100ps-1800ext/compiled/defect-cut-off-{defect_cut_off}A/reshaped")
compiled_data.mkdir(exist_ok=True)

lipid_compositions = {
    "DOPC": [0, 0],
    "DPPC": [0, 1],
    "DOPC-DOPA": [1, 0],
    "DPPC-DPPA": [1, 1],
    "DOPC-DOPS": [2, 0],
    "DPPC-DPPS": [2, 1]
}


# lipid_compositions = {
#     "DOPC": [0, 0],
#     "DOPC-DOPA": [1, 0],
#     "DOPC-DOPS": [2, 1]
# }


# lipid_compositions = {
#     "DPPC": [0, 0],
#     "DPPC-DPPA": [1, 0],
#     "DPPC-DPPS": [2, 1]
# }

# lipid_compositions = {
#     "DOPC-DOPS": [0, 0], 
#     "DPPC-DPPS": [0, 1]
# }



lipid_comp = list(lipid_compositions.keys())



############################################################################
#bar plot avg defect size over frames for a restricted region of membrane ###
#############################################################################

# def bar_graph_system_peak(data,restriction,configuration):
#     systems = list(data.keys())
#     avg_defects_coverage = [data[sys][0] for sys in systems]
#     std_defect_coverage = [data[sys][1] for sys in systems]

#     #set labels for system
#     x = np.arange(len(systems))  

#     width = 0.35 
    
#     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
#         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=3)
#         #bar1 = ax.bar(x - width/2, avg_peak_defects, width, yerr=std_peak_defects, label='Peak Defects', capsize=5)
#         bar2 = ax.bar(x + width/2,  avg_defects_coverage, width ,yerr = std_defect_coverage, label='Total Defects', capsize=5)

#         ax.set_xlabel('Systems')
#         ax.set_ylabel('Average Surface Area Covered by Defects (%)')
#         ax.set_xticks(x)
#         ax.set_ylim(0,10.5)

#         ax.set_xticklabels(lipid_comp, rotation=45, ha='right')

#         # Add value annotations on top of bars
#         for bar_group in [bar2]:
#             for bar in bar_group:
#                 height = bar.get_height()
#                 # ax.annotate(f'{height:.3f}',
#                 #             xy=(bar.get_x() + bar.get_width() / 2, height),
#                 #             xytext=(0, 3), 
#                 #             textcoords="offset points",
#                 #             ha='center', va='bottom')


#         fig.tight_layout()
#         plt.savefig(f"defect-coverage-bar-plt-avg-frames-{configuration}-restriction{restriction}.png")


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
    system_spacing = 4   # separates systems
    config_spacing = 1  # keeps configs tight
    
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.2, subheight=1.5)
        

        #plot each configuration
        for i, config in enumerate(configurations):
            avg_values = [data[sys][i][0] /100 for sys in systems] #/100 to covert to nm^2
            std_values = [data[sys][i][1] /100 for sys in systems] #/100 to covert to nm^2
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

        plt.savefig(f"{compiled_data_reshaped }/reshaped-{config}-peak{restriction}-{lipid_type}-cutoff{cut_off}-defectsize.png")
        plt.savefig(f"{compiled_data_reshaped }/reshaped-{config}-peak{restriction}-{lipid_type}-cutoff{cut_off}-defectsize.pdf")

#######################################################
##execution code for bar graphs for avg defect size ###
#######################################################
restriction = 40
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
            "avg_defect_size": avg,
            "std": std,
            "num_frames": num_frames
        })
df = pd.DataFrame(rows)
df.to_csv(f"{compiled_data_reshaped}/reshaped_{lipid_type}_defect_data_avg_size{restriction}-cutoff{cut_off}.csv", index = False)

#######################################################
##execution code for bar graphs for defect coverage ###
#######################################################
restriction = 40
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
            "avg_defect_coverage": avg,
            "std": std,
            "num_frames": num_frames
        })
df = pd.DataFrame(rows)
df.to_csv(f"{compiled_data_reshaped }/reshaped_{lipid_type}_defect_data_avg_coverage{restriction}_cutoff{cut_off}.csv", index = False)



# # # ################################################
# # # #scatter plot for defect coverage over frames###
# # # ################################################
# # goal: to see if the surface defect coverage converges
# #scatter plot for defect coverage over frames
def defect_coverage_simulation(defect_restricted_per_frame_system, avg_defect_DPPC, avg_defect_DOPC):
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):

        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=5, subheight=3)

        configuration = "flat"
        restriction = 40

        colors_presets = [ "#f8c362ff",  "#5db7fc" ]


        for i, composition in enumerate(lipid_compositions):

            per_frame_array = \
                defect_restricted_per_frame_system[composition][configuration][restriction]

            frames = np.arange(len(per_frame_array))

            ax.scatter(
                frames,
                per_frame_array,
                s=3,
                color=colors_presets[i],
                label=composition
            )
        ax.axhline(y=avg_defect_DOPC, color = "black") #DOPC 
        ax.axhline(y=avg_defect_DPPC, color = "red") #DPPC
        ax.set_ylim(0,22)
        #ax.set_xlim(400,)
        ax.set_xlabel("Frame")
        ax.set_ylabel("Percent defect coverage")
        ax.set_title(f"Defect coverage per frame ({configuration}, restriction={restriction} Å)")

        ax.legend(frameon=False)

        # plt.tight_layout()

        # plt.savefig(
        #     f"newfig.png",
        #     dpi=40
        # )

        # plt.savefig(
        #     f"frame-defect-coverage-{configuration}-entire-restriction{restriction}-flat-1800ext.png",
        #     dpi=40
        # )


# # ###########################################################
# # ###calculate the size of the defects in a defined region### 
# # ###########################################################
# restriction = 40
# system_defect_sizes = {}
# for composition, position in lipid_compositions.items():
#     for configuration in configurations: 
#         #locate the appropriate data folder
#         if configuration ==  "flat": 
#             time = "100ps-1800ext"
#             analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
#             analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
#             folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
#         if configuration == "strain2":
#             time = "100ps-1800ext"
#             analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
#             analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"
#     defect_counts = np.load(folder_system/"output-defects-upper.npy")
#     quad_centroids_np = np.load(folder_system/"output-centroids.npy")

#     all_defect_sizes = [] #a list that states each defect size
#     n_frames = len(quad_centroids_np)

#     for frame in np.arange(0, n_frames):
#         defects_counts_frame =  defect_counts[frame]
#         centroids_frame =  quad_centroids_np[frame, :, 0:1]
#         #restrain to region of interest 
#         centered_centroids = centroids_frame  - np.mean(centroids_frame)

#         centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()

#         restricted_region_grids_with_defects =  np.where((centroid_X_mask == 1) & ( defects_counts_frame  < 1),1,0)

#         #need to reshape numpy before apply the scipy 
#         label_all, nfeat_all = scipy_ndimage.label(restricted_region_grids_with_defects)

#         for i in range(1, nfeat_all + 1):
#             defect_size = np.count_nonzero(label_all == i)
#             all_defect_sizes.append(defect_size)

#     system_defect_sizes[composition] = all_defect_sizes 


# # # # #############################################################
# # # # #####bar graph to show the distributions of defect sizes#####
# # # # #############################################################
# # each bin is 2A, to help with visualization
# restriction = 40
# system_defect_sizes = {}
# system_defect_sizes_log_transformed_count = {}
# for composition, position in lipid_compositions.items():
#     strain_defect_sizes_log_transformed_count = {}
#     for configuration in configurations:
#         if configuration == "flat":
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
#         if configuration == "strain2":
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"
        
#         defect_counts = np.load(folder_system/"output-defects-upper.npy")
#         quad_centroids_np = np.load(folder_system/"output-centroids.npy")
#         n_frames = len(quad_centroids_np)

#         defect_sizes = []
#         defect_size_counts = []
#         defect_size_counts_log_transformed = []
#         bins_below50 = np.arange(0, 50, 2) #each bin is 2A
#         bins_above50 = np.arange(50, 451, 400)
#         bins = np.concatenate([bins_below50, bins_above50])
#         for frame in range(0, n_frames):
#             defects_counts_frame =  defect_counts[frame]

#             #all defect mask
#             masked_all_defects_num = np.where(
#                 defects_counts_frame < 1,
#                 1,
#                 0
#             )

#             centroids_frame =  quad_centroids_np[frame, :, 0:1]

#             #restrain to region of interest 
#             centered_centroids = centroids_frame  - np.mean(centroids_frame)
#             centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()
#             restricted_region_grids_with_defects =  np.where((centroid_X_mask == 1) & (masked_all_defects_num == 1),1,0)

#             #quantifying the size of the defect 
#             #need to reshape numpy before applying the scipy
#             label_all, nfeat_all = scipy_ndimage.label( restricted_region_grids_with_defects)
#             for i in range(1, nfeat_all + 1):
#                 defect_size = np.count_nonzero(label_all == i)
#                 defect_sizes.append(defect_size)
        
#         print(composition, "max defect size:",np.max(defect_sizes) )
#         #bin the data
#         data1, _ = np.histogram(
#             defect_sizes, bins = bins
#         )
#         bin_centers = 0.5 * (bins[1:] + bins[:-1]) 

#         #log transform
#         log_transform = np.log(data1)

#         bin_widths = np.diff(bins)

#         #saves the log_transformed data per configuration of the same lipid composition type
#         strain_defect_sizes_log_transformed_count[configuration]=log_transform

#     system_defect_sizes_log_transformed_count[composition] = strain_defect_sizes_log_transformed_count



# # for configuration in configurations:
# #     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
# #         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.40, subheight=1.40)
# #         #colors_presets = [ "#5db7fc","#f8c362ff"]
# #         #colors_presets = [ "#d2826cff","#aaa9a9ff"]
# #         #i = 0
# #         alpha = 1
# #         for composition in lipid_comp:
# #             #color = colors_presets[i]
# #             ax.bar(bin_centers, system_defect_sizes_log_transformed_count[composition][configuration], alpha = alpha, width = bin_widths, label = composition)
# #             ax.set_ylabel("log(defect-counts)")
# #             ax.xaxis.set_major_locator(MultipleLocator(10))
# #             ax.set_xlabel("Defect Area (Å²)")
# #             alpha = alpha 
# #             ax.set_xlim(0,55)
# #             ax.set_title(configuration)
# #             #i += 1
# #             alpha -= .40
# #     # plt.legend()
# #     # plt.savefig(f"{compiled_data}/logdefect-barplot-fig-peak-{lipid_type}-{configuration}-ext1800-cutoff.pdf", dpi = 40)
# #     # plt.savefig(f"{compiled_data}/logdefect-barplot-fig-peak-{lipid_type}-{configuration}-ext1800-cutoff.png", dpi = 40)
# #     # plt.close()

# # ###########################################################
# # ###calculate the size of the defects in a defined region### 
# # ###########################################################
# restriction = 40
# system_defect_sizes = {}
# for composition, position in lipid_compositions.items():
#     for configuration in configurations: 
#         #locate the appropriate data folder
#         if configuration ==  "flat": 
#             time = "100ps-1800ext"
#             analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
#             analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
#             folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
#         if configuration == "strain2":
#             time = "100ps-1800ext"
#             analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
#             analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"
#     defect_counts = np.load(folder_system/"output-defects-upper.npy")
#     quad_centroids_np = np.load(folder_system/"output-centroids.npy")

#     all_defect_sizes = [] #a list that states each defect size
#     n_frames = len(quad_centroids_np)

#     for frame in np.arange(0, n_frames):
#         defects_counts_frame =  defect_counts[frame]
#         centroids_frame =  quad_centroids_np[frame, :, 0:1]
#         #restrain to region of interest 
#         centered_centroids = centroids_frame  - np.mean(centroids_frame)

#         centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()

#         restricted_region_grids_with_defects =  np.where((centroid_X_mask == 1) & ( defects_counts_frame  < 1),1,0)

#         #need to reshape before applying the scipy 
#         label_all, nfeat_all = scipy_ndimage.label( restricted_region_grids_with_defects)

#         for i in range(1, nfeat_all + 1):
#             defect_size = np.count_nonzero(label_all == i)
#             all_defect_sizes.append(defect_size)

#     system_defect_sizes[composition] = all_defect_sizes 


# # # # #############################################################
# # # # #####line-graph to show the distributions of defect sizes#####
# # # # #############################################################
# # each bin is 1A 
# restriction = 40
# system_defect_sizes_log_transformed_count = {}
# for composition, position in lipid_compositions.items():
#     strain_defect_sizes_log_transformed_count = {}
#     for configuration in configurations:
#         if configuration == "flat":
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
#         if configuration == "strain2":
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"
        
#         defect_counts = np.load(folder_system/"output-defects-upper.npy")
#         quad_centroids_np = np.load(folder_system/"output-centroids.npy")
#         n_frames = len(quad_centroids_np)
#         shape = np.load(folder_system/"output-shape.npy")

#         defect_sizes = []
#         defect_size_counts = []
#         defect_size_counts_log_transformed = []
#         # bins_below50 = np.arange(0, 50 1) #each bin is 2A
#         # bins_above50 = np.arange(50, 451, 400)
#         # bins_below100 = np.arange(0, 40, 1) #each bin is 2A
#         # bins_above100 = np.arange(40, 700, 400)
#         bins = np.arange(0,4000,1)
#         #bins = np.concatenate([bins_below100, bins_above100])
#         for frame in range(0, n_frames):
#             defects_counts_frame =  defect_counts[frame]
#             defects_counts_frame = defects_counts_frame.reshape(shape)

#             #all defect mask
#             masked_all_defects_num = np.where(
#                 defects_counts_frame < 1,
#                 1,
#                 0
#             )

#             centroids_frame =  quad_centroids_np[frame, :, 0:1]
            
#             #restrain to region of interest 
#             centered_centroids = centroids_frame  - np.mean(centroids_frame)
#             centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()
#             centroid_mask_grid = centroid_X_mask.reshape(shape)

#             restricted_region_grids_with_defects =  np.where((centroid_mask_grid == 1) & (masked_all_defects_num == 1),1,0)

#             #quantifying the size of the defect 
#             #need to reshape 
#             label_all, nfeat_all = scipy_ndimage.label( restricted_region_grids_with_defects)
#             for i in range(1, nfeat_all + 1):
#                 defect_size = np.count_nonzero(label_all == i)
#                 if defect_size >=15:
#                     defect_sizes.append(defect_size)
        
#         print(composition, "max defect size:",np.max(defect_sizes) )
#         #bin the data
#         data1, _ = np.histogram(
#             defect_sizes, bins = bins
#         )
#         bin_centers = bins[:-1]


#         #log transform
#         log_transform = np.log(data1)



#         #saves the log_transformed data per configuration of the same lipid composition type
#         strain_defect_sizes_log_transformed_count[configuration]=log_transform

#     system_defect_sizes_log_transformed_count[composition] = strain_defect_sizes_log_transformed_count


# # for configuration in configurations:
# #     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
# #         #fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.40, subheight=1.40)
# #         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4.40, subheight=4.40)
# #         colors_presets = [ "#ff6c0A","#8ece33", "#9A36D8"]
# #         i = 0
# #         alpha = 1
# #         for composition in lipid_comp:
# #             color = colors_presets[i]
# #             ax.scatter(bin_centers, system_defect_sizes_log_transformed_count[composition][configuration], alpha = alpha, label = composition, s =1)
# #             ax.set_ylabel("log(defect-counts)")
# #             ax.xaxis.set_major_locator(MultipleLocator(400))
# #             ax.set_xlabel("Defect Area (Å²)")
# #             alpha = alpha 
# #             #ax.set_xlim(15,55)
# #             #ax.set_xlim(140,1000)
# #             ax.set_ylim(0,10)
# #             ax.set_title(configuration)
# #             i += 1
# #     plt.legend()
# #     plt.savefig(f"{compiled_data}/reshape-logdefect-linegraph-fig-peak-{lipid_type}-{configuration}-retriction{restriction}-ext1800-cutoff-cutoff{cut_off}.pdf", dpi = 40)
# #     plt.savefig(f"{compiled_data}/reshape-logdefect-linegraph-fig-peak-{lipid_type}-{configuration}-retriction{restriction}-ext1800-cutoff-cutoff{cut_off}.png", dpi = 40)
# #     plt.close()




# #####################################################################
# ###plot for the probability of finding a defect of a certain size### 
# #####################################################################
# major_probs = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]) 
# restriction = 40
# system_defect_sizes_log_transformed_prob = {}
# for composition, position in lipid_compositions.items():
#     strain_defect_sizes_log_transformed_prob = {}
#     for configuration in configurations:
#         # create a new figure for each system

#         if configuration == "flat":
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
#         if configuration == "strain2":
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"

#         defect_counts = np.load(folder_system/"output-defects-upper.npy")
#         quad_centroids_np = np.load(folder_system/"output-centroids.npy")
#         n_frames = len(quad_centroids_np)
#         shape = np.load(folder_system/"output-shape.npy")

#         defect_sizes = []
#         bins = np.arange(0,4000,1)
#         # bins_below40 = np.arange(0, 40, 1) #each bin is 1A
#         # bins_above40 = np.arange(40, 4000, 1 )
#         # bins = np.concatenate([bins_below40, bins_above40])

#         for frame in range(n_frames):

#             defects_counts_frame = defect_counts[frame].reshape(shape)

#             masked_all_defects_num = np.where(defects_counts_frame < 1, 1, 0)

#             centroids_frame = quad_centroids_np[frame, :, 0:1]

#             centered_centroids = centroids_frame - np.mean(centroids_frame)
#             centroid_X_mask = np.where(
#                 (centered_centroids >= -restriction) & 
#                 (centered_centroids <= restriction), 1, 0
#             ).squeeze()

#             centroid_mask_grid = centroid_X_mask.reshape(shape)

#             restricted_region_grids_with_defects = np.where(
#                 (centroid_mask_grid == 1) & 
#                 (masked_all_defects_num == 1), 1, 0
#             )

#             label_all, nfeat_all = scipy_ndimage.label(restricted_region_grids_with_defects)

#             for i in range(1, nfeat_all + 1):
#                 defect_size = np.count_nonzero(label_all == i)
#                 if defect_size >= 15:
#                     defect_sizes.append(defect_size)

#             # histogram
#             counts, _ = np.histogram(defect_sizes, bins=bins)
#             probabilities = counts / counts.sum()
#             bin_centers = 0.5 * (bins[1:] + bins[:-1])

#             mask = probabilities > 0
#             areas_nonzero = bin_centers[mask]
#             probabilities_nonzero = probabilities[mask]
#             log_probabilities = np.log(probabilities_nonzero)

#             # apply PackMem-like fit limits
#             masks_limits = (
#                 (areas_nonzero >= lower_limit_size_defect) &
#                 (probabilities_nonzero >= 1e-4)
#             )


#             areas_nonzero = areas_nonzero[masks_limits]
#             probabilities_nonzero = probabilities_nonzero[masks_limits]
#             log_probabilities = np.log(probabilities_nonzero)


#             #saves the log probability data per configuration 
#             strain_defect_sizes_log_transformed_prob[configuration]= [areas_nonzero,log_probabilities]
        
#         system_defect_sizes_log_transformed_prob[composition] = strain_defect_sizes_log_transformed_prob



# for configuration in configurations:
#     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
#         #fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.3, subheight=1.4)
#         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=2.2, subheight=1.5) #large graphs
#         colors_presets = [ "#ff6c0A","#8ece33", "#9A36D8"]
#         i = 0
#         for composition in lipid_comp:
#             color = colors_presets[i]
#             areas_nonzero = system_defect_sizes_log_transformed_prob[composition][configuration][0]
#             log_probabilities = system_defect_sizes_log_transformed_prob[composition][configuration][1]
#             # scatter
#             ax.scatter(areas_nonzero ,log_probabilities , s=1, label = composition, color = color)
#             i += 1

#         # # linear regression
#         # linear_regression = linregress(areas_nonzero, log_probabilities)

#         # slope = linear_regression.slope
#         # intercept = linear_regression.intercept
#         # r2_value = linear_regression.rvalue**2
#         # packdef_constant = abs(1 / slope)

#         # xfit = np.linspace(min(areas_nonzero), max(areas_nonzero), 100)
#         # yfit = slope * xfit + intercept

#         # ax.plot(xfit, yfit, color="black",
#         #         label=f"y={slope:.2f}x+{intercept:.2f}\n$R^2$={r2_value:.2f}\n defect-constant={packdef_constant:.2f}")

#         # axis formatting
#         ax.set_yticks(np.log(major_probs))
#         ax.set_yticklabels([
#             r"$10^{-1}$", r"$10^{-2}$", r"$10^{-3}$",
#             r"$10^{-4}$", r"$10^{-5}$", r"$10^{-6}$"
#         ])

#         ax.set_ylabel("Probability")
#         ax.set_xlabel("Defect Area (Å²)")
#         ax.set_ylim(np.log(5.e-5),-1)
#         #ax.set_xlim(0,1100) #strain limit
#         #ax.set_xlim(0,40) #flat limit 
#         ax.set_xlim(0,400) #flat limit 


#         ax.set_title(f"{lipid_type} — {configuration}")
#         plt.legend()

#         plt.tight_layout()

#         plt.savefig(
#             f"{compiled_data_reshaped}/reshaped-{lipid_type}-{configuration}-defect-size-probability.png"
#         )
#         plt.savefig(
#             f"{compiled_data_reshaped}/reshaped-{lipid_type}-{configuration}-defect-size-probability.pdf"
#         )
#         plt.close()
    