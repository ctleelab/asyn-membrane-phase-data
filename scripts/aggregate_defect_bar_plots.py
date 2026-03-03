

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
from scipy import ndimage as scipy_ndimage
from math import log
from collections import Counter


configurations = ["flat"]
defect_cut_off = 10
#time = "100ps-1200ext"
defect_type = ["deep", "shallow", "all"]
lower_limit_size_defect = 15
upper_limit_size_defect = 60
restrictions = [200]

# lipid_compositions = {
#     "DOPC": [0, 0],
#     "DPPC": [0, 1],
#     "DOPC-DOPA": [1, 0],
#     "DPPC-DPPA": [1, 1],
#     "DOPC-DOPS": [2, 0],
#     "DPPC-DPPS": [2, 1]
# }

lipid_compositions = {
    "DOPC": [0, 0],
    "DPPC": [0, 1]}



lipid_comp = list(lipid_compositions.keys())



# #############################################################################
# # #bar plot avg defect size over frames for a restricted region of membrane ###
# # #############################################################################

# # def bar_graph_system_peak(data,restriction,configuration):
# #     systems = list(data.keys())
# #     avg_defects_coverage = [data[sys][0] for sys in systems]
# #     std_defect_coverage = [data[sys][1] for sys in systems]

# #     #set labels for system
# #     x = np.arange(len(systems))  

# #     width = 0.35 
    
# #     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
# #         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=3)
# #         #bar1 = ax.bar(x - width/2, avg_peak_defects, width, yerr=std_peak_defects, label='Peak Defects', capsize=5)
# #         bar2 = ax.bar(x + width/2,  avg_defects_coverage, width ,yerr = std_defect_coverage, label='Total Defects', capsize=5)

# #         ax.set_xlabel('Systems')
# #         ax.set_ylabel('Average Surface Area Covered by Defects (%)')
# #         ax.set_xticks(x)
# #         ax.set_ylim(0,10.5)

# #         ax.set_xticklabels(lipid_comp, rotation=45, ha='right')

# #         # Add value annotations on top of bars
# #         for bar_group in [bar2]:
# #             for bar in bar_group:
# #                 height = bar.get_height()
# #                 # ax.annotate(f'{height:.3f}',
# #                 #             xy=(bar.get_x() + bar.get_width() / 2, height),
# #                 #             xytext=(0, 3), 
# #                 #             textcoords="offset points",
# #                 #             ha='center', va='bottom')


# #         fig.tight_layout()
# #         plt.savefig(f"defect-coverage-bar-plt-avg-frames-{configuration}-restriction{restriction}.png")


def bar_graph_system_peak_grouped(data,restriction,configurations):
    systems = list(data.keys())
    n_configs = len(configurations)

    #set labels for system
    x = np.arange(len(systems))  

    width = 0.35 
    system_spacing = 4   # separates systems
    config_spacing = 1  # keeps configs tight
    
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=2)
        

        #plot each configuration
        for i, config in enumerate(configurations):
            avg_values = [data[sys][i][0] for sys in systems]
            std_values = [data[sys][i][1] for sys in systems]
            group_spacing = 1.2 
            offsets = (i - (n_configs - 1) / 2) * width * config_spacing
            print(systems)
            print(f"{config}-avg-defect-coverage:", avg_values)

            ax.bar(x + offsets,  avg_values, width ,yerr = std_values, label=config, capsize=5)
        ax.legend()
        ax.set_ylabel('Surface Covered by Defects (%)')
        ax.set_xticks(x)
        ax.set_ylim(0,20)

        ax.set_xticklabels(lipid_comp, rotation=45, ha='right')

        plt.savefig(f"all-systems-peak{restriction}.png")


def total_num_defect_restricted_X(defect_data,shape,centroids,restriction,configuration):
    
    #center of grids at 0 
    centered_centroids = centroids - np.mean(centroids)

    centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()

    restricted_region_grids_with_defects =  np.where((centroid_X_mask == 1) & (defect_data == 1),1,0)

    total_grids = len(defect_data)


    defect_grids = np.sum(restricted_region_grids_with_defects)



    percent_defect_coverage = (defect_grids/total_grids)*100

    #return total_num_defects_restricted_X, total_grids_restricted_X
    return defect_grids,  total_grids,  percent_defect_coverage 


####Surface covered by defects ####
def surface_coverage_defects(restrictions): 
    defect_restricted_avg_system = {} 
    defect_restricted_per_frame_system = {}
    frames_avg_defect_coverage_per_system = {}
    for composition in lipid_comp:
        configurations_types = {}
        configurations_types_per_frame = {}
        configurations_types_frames_to_avg_defect = {}
        for configuration in configurations: 
            if configuration ==  "flat": 
                time = "100ps-1200ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext"
            if configuration == "strain2":
                time = "100ps-1800ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800"
            print("configuration:", configuration)
            defect_counts = np.load(folder_system/"output-defects-upper.npy")
            quad_centroids_np = np.load(folder_system/"output-centroids.npy")
            shape = np.load(folder_system/"output-shape.npy")
            X = np.load(folder_system/"output-Xshape.npy")
            Y = np.load(folder_system/"output-Yshape.npy")
            print("shape-defects:", defect_counts.shape)
            n_frames = len(quad_centroids_np)

            percent_cut_defects = {}
            percent_cut_defects_per_frame  = {}
            for restriction in restrictions: 

                total_defects_restricted_region_per_frame = np.empty(n_frames)
                total_grids_restricted_region_per_frame = np.empty(n_frames)
                percentage_defect_coverage_per_frame = np.empty(n_frames)
                list_frames_close_avg_defect = []


                for frame in np.arange(0, n_frames):


                    defect_counts_frame = defect_counts[frame]
                    centroids_frame =  quad_centroids_np[frame, :, 0:1]

                    #all defect 
                    masked_all_defects_num = np.where(
                        defect_counts_frame < 1,
                        1,
                        0
                    )

                    total_per_frame, total_grids, percent_defect_coverage = total_num_defect_restricted_X(masked_all_defects_num,shape,centroids_frame,restriction,configuration)

                    total_defects_restricted_region_per_frame[frame] = total_per_frame

                    total_grids_restricted_region_per_frame [frame] = total_grids

                    percentage_defect_coverage_per_frame[frame] = percent_defect_coverage

                average_defect_coverage = np.mean(percentage_defect_coverage_per_frame)
                std_defect_coverage = np.std(percentage_defect_coverage_per_frame)

                percent_cut_defects[restriction]  = [average_defect_coverage,std_defect_coverage]

                #storing per frame array, to look at convergence
                percent_cut_defects_per_frame[restriction] = percentage_defect_coverage_per_frame

            configurations_types[configuration] = percent_cut_defects
            configurations_types_per_frame[configuration] = percent_cut_defects_per_frame
            configurations_types_frames_to_avg_defect[configuration] = list_frames_close_avg_defect
        defect_restricted_avg_system[composition]= configurations_types
        defect_restricted_per_frame_system[composition] = configurations_types_per_frame
        frames_avg_defect_coverage_per_system[composition] = configurations_types_frames_to_avg_defect

    return defect_restricted_avg_system
defect_restricted_avg_system = surface_coverage_defects(restrictions)

def frame_extraction_avg_surface_defect_coverage(defect_restricted_avg_system):
    '''
    Compiles a list of frames that is +/- .01 away from the avg defect coverage over the entire membrane

    Args:

    Returns:
    max_frames_avg_defect_coverage_per_system (dict): A dictionary that contains a singular frame that represents the avg coverage of defects across the entire membrane for each lipid composition.
    '''
    frames_avg_defect_coverage_per_system = {} 
    max_frames_avg_defect_coverage_per_system = {}
    for composition in lipid_comp:
        configurations_types_frames_to_avg_defect = {}
        list_frames_close_avg_defect = []
        for configuration in configurations:
            if configuration ==  "flat": 
                time = "100ps-1200ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path  / f"{composition}-{configuration}-production-stripped-ext"
            if configuration == "strain2":
                time = "100ps-1800ext"
                analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}"
                analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
                folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped-centered-ext1800"

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

                    _, _, percent_defect_coverage = total_num_defect_restricted_X(masked_all_defects_num,shape,centroids_frame,restriction,configuration)

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
        configurations_types_frames_to_avg_defect[configuration] = list_frames_close_avg_defect
        frames_avg_defect_coverage_per_system[composition] = configurations_types_frames_to_avg_defect
        max_frames_avg_defect_coverage_per_system[composition] = max(frames_avg_defect_coverage_per_system[composition][configuration]) 
    return max_frames_avg_defect_coverage_per_system

max_frames_avg_defect_coverage_per_system = frame_extraction_avg_surface_defect_coverage(defect_restricted_avg_system)
print(max_frames_avg_defect_coverage_per_system)

restriction = 200
data_grouped = {}

for system in lipid_comp:
    data_grouped[system] = []
    for config in configurations:
        avg, std = defect_restricted_avg_system[system][config][restriction]
        data_grouped[system].append([avg, std])

    #bar_graph_system_peak(temp_dict, restriction)

bar_graph_system_peak_grouped(data_grouped, restriction, configurations)



# # # ################################################
# # # #scatter plot for defect coverage over frames###
# # # ################################################
# #scatter plot for defect coverage over frames
# with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):

#     fig, ax = ph.fixed_size_subplots(1, 1, subwidth=5, subheight=3)

#     configuration = "flat"
#     restriction = 200

#     colors_presets = [ "#f8c362ff",  "#5db7fc" ]


#     for i, composition in enumerate(lipid_compositions):

#         per_frame_array = \
#             defect_restricted_per_frame_system[composition][configuration][restriction]

#         frames = np.arange(len(per_frame_array))

#         ax.scatter(
#             frames,
#             per_frame_array,
#             s=3,
#             color=colors_presets[i],
#             label=composition
#         )
#     ax.axhline(y=avg_defect_DOPC, color = "black") #DOPC 
#     ax.axhline(y=avg_defect_DPPC, color = "red") #DPPC
#     ax.set_ylim(0,22)
#     #ax.set_xlim(3000,)
#     ax.set_xlabel("Frame")
#     ax.set_ylabel("Percent defect coverage")
#     ax.set_title(f"Defect coverage per frame ({configuration}, restriction={restriction} Å)")

#     ax.legend(frameon=False)

#     plt.tight_layout()

#     plt.savefig(
#         f"newfig.png",
#         dpi=300
#     )

#     plt.savefig(
#         f"frame-defect-coverage-{configuration}-entire-restriction{restriction}-PS-flat-1200ext.png",
#         dpi=300
#     )



# # # ###########################################################
# # # ###calculate the size of the defects in a defined region### 
# # # ###########################################################

# # restriction = 200
# # system_defect_sizes = {}
# # for composition, position in lipid_compositions.items():
# #     folder_system = analysis_defect_path / f"{composition}-{configurations}-production-stripped_centered"
# #     defect_counts = np.load(folder_system/"output-defects-upper.npy")
# #     quad_centroids_np = np.load(folder_system/"output-centroids.npy")

# #     all_defect_sizes = [] #a list that states each defect size
# #     n_frames = len(quad_centroids_np)


# #     for frame in np.arange(0, n_frames):
# #         defects_counts_frame =  defect_counts[frame]
# #         centroids_frame =  quad_centroids_np[frame, :, 0:1]
# #         #restrain to region of interest 
# #         centered_centroids = centroids_frame  - np.mean(centroids_frame)

# #         centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()

# #         restricted_region_grids_with_defects =  np.where((centroid_X_mask == 1) & ( defects_counts_frame  < 1),1,0)


# #         label_all, nfeat_all = scipy_ndimage.label( restricted_region_grids_with_defects)

# #         for i in range(1, nfeat_all + 1):
# #             defect_size = np.count_nonzero(label_all == i)
# #             all_defect_sizes.append(defect_size)

# #     system_defect_sizes[composition] = all_defect_sizes 


# # ##########################################
# # ###bar plot avg defect size over frames### 
# # ##########################################
# # def bar_graph_system(data):
# #     systems = list(data.keys())
# #     list_of_defect_sizes  = [data[sys] for sys in systems]
# #     avg_defects_size = [np.mean(sizes) for sizes in list_of_defect_sizes]
# #     std_defect_size = [np.std(sizes) for sizes in list_of_defect_sizes]



# #     #set labels for system
# #     x = np.arange(len(systems))  

# #     width = 0.35 
    
# #     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
# #         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=3)
# #         #bar1 = ax.bar(x - width/2, avg_peak_defects, width, yerr=std_peak_defects, label='Peak Defects', capsize=5)
# #         bar2 = ax.bar(x + width/2,  avg_defects_size, width, yerr = std_defect_size, label='Total Defects', capsize=5)

# #         ax.set_xlabel('Systems')
# #         ax.set_ylabel('Average defect size A^2')
# #         ax.set_xticks(x)

# #         ax.set_xticklabels(lipid_comp, rotation=45, ha='right')
# #         ax.legend()

# #         # Add value annotations on top of bars
# #         for bar_group in [bar2]:
# #             for bar in bar_group:
# #                 height = bar.get_height()
# #                 ax.annotate(f'{height:.1f}',
# #                             xy=(bar.get_x() + bar.get_width() / 2, height),
# #                             xytext=(0, 3), 
# #                             textcoords="offset points",
# #                             ha='center', va='bottom')


# #         fig.tight_layout()
# #         plt.savefig(f"defect-size-bar-plt-avg-frames-{configurations}.png")
# #         plt.show()

# # # bar_graph_system(system_defect_sizes )



# # # ###################################
# # # ###distrubution of defect sizes ### 
# # # ###################################
# # # #area under curve is 1, so shape of defect histograms are comparable

# # system_defect_sizes = {}
# # with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
# #     fig,ax = ph.fixed_size_subplots(1, 1, subwidth=1, subheight=1)
# #     for composition, position in lipid_compositions.items():
# #         folder_system = analysis_defect_path / f"{composition}-{configurations}-production-stripped_centered"
# #         defect_counts = np.load(folder_system/"output-defects-upper.npy")

# #         all_defect_sizes = [] #a list that states has each defect size 
# #         bins = np.arange(0, 400, 1) #define bins so that they are equal between systems

# #         for frame_defects in defect_counts:
# #             masked_all_defects_num = np.where((frame_defects < 1), 1, 0)

# #             label_all, nfeat_all = scipy_ndimage.label(masked_all_defects_num)

# #             for i in range(1, nfeat_all + 1):
# #                 defect_size = np.count_nonzero(label_all == i)
# #                 all_defect_sizes.append(defect_size)
# #         system_defect_sizes[composition] = all_defect_sizes 


# #         pdf, bins = np.histogram(all_defect_sizes, bins=bins, density=True)
# #         bin_centers = 0.5 * (bins[1:] + bins[:-1])


# #         ax.bar(bin_centers, pdf, label = composition)

# # ax.set_ylabel("?")
# # ax.set_xlim(0,55)
# # ax.set_xlabel("Defect Area (A^2)")
# # ax.legend(loc = "upper right", fontsize = "small")
# # plt.tight_layout()
# # plt.savefig(f"all-probability-{composition}-defect-size-{configurations}")
# # plt.close()



# # # # ############################################
# # # # ###calculate the maximum size of a defect### 
# # # # ############################################

# # largest_defect = 0
# # restriction = 200
# # system_largest_defect_sizes_max = {}
# # system_largest_defect_sizes_max_per_frame = {}
# # for composition, position in lipid_compositions.items():
# #     largest_defect = 0
# #     folder_system = analysis_defect_path / f"{composition}-{configurations}-production-stripped_centered"
# #     defect_counts = np.load(folder_system/"output-defects-upper.npy")
# #     quad_centroids_np = np.load(folder_system/"output-centroids.npy")
# #     n_frames = len(quad_centroids_np)

# #     max_defect_sizes = np.empty(n_frames) 
# #     n_frames = len(quad_centroids_np)


# #     for frame in range(n_frames):
# #         defects_counts_frame =  defect_counts[frame]

# #         #all defect 
# #         masked_all_defects_num = np.where(
# #             defects_counts_frame < 1,
# #             1,
# #             0
# #         )

# #         centroids_frame =  quad_centroids_np[frame, :, 0:1]
# #         #restrain to region of interest 
# #         centered_centroids = centroids_frame  - np.mean(centroids_frame)

# #         centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()

# #         restricted_region_grids_with_defects =  np.where((centroid_X_mask == 1) & (masked_all_defects_num == 1),1,0)


# #         label_all, nfeat_all = scipy_ndimage.label( restricted_region_grids_with_defects)

        
# #         if nfeat_all == 0:
# #             max_defect_sizes[frame] = 0
# #         else:
# #             defect_sizes = []
# #             for i in range(1, nfeat_all + 1):
# #                 defect_sizes.append(np.count_nonzero(label_all == i))

# #             max_defect_sizes[frame] = np.max(defect_sizes)

# #     system_largest_defect_sizes_max_per_frame[composition] = max_defect_sizes
                    

# #     #system_largest_defect_sizes_max[composition] = largest_defect 
# # print(system_largest_defect_sizes_max_per_frame)



# # #scatter plot for defect_restricted_system_per_frame
# # plt.figure(figsize=(10, 6))

# # with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
# #     fig, ax = ph.fixed_size_subplots(2, 1, subwidth=3, subheight=3)
# #     ax = np.array(ax).flatten() 
# #     for idx, (composition, position) in enumerate(lipid_compositions.items()):

# #         # Extract the per-frame defect coverage array
# #         per_frame_array = system_largest_defect_sizes_max_per_frame[composition]
# #         print(per_frame_array.shape)

# #         frames = np.arange(0,6001)
# #         # Scatter plot
# #         ax[idx].scatter(frames, per_frame_array, s =3, label=f"{composition}-{configurations}", alpha=0.6)
# #         ax[idx].set_ylim(200.5,100)
# #         ax[idx].set_xlabel("Frame")
# #         ax[idx].set_ylabel("largest defect A^2")
# #         ax[idx].set_title(composition)

# # # plt.xlabel("Frame")
# # # plt.ylim(200.5,)
# # # plt.ylabel("largest defect A^2")
# # # plt.title(f"Defect coverage per frame (restriction={restriction})")
# # # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
# # plt.tight_layout()
# # plt.savefig(f"frame-largest-defect-{configurations}-restriction{restriction}-PC.png")










# # # # #############################################
# # # # #####bar graph to show the distributions#####
# # # # #############################################

# restriction = 200
# system_defect_sizes = {}
# system_defect_sizes_log_transformed_count = {}
# for composition, position in lipid_compositions.items():
#     folder_system = analysis_defect_path / f"{composition}-{configurations}-production-stripped-centered-ext"
#     defect_counts = np.load(folder_system/"output-defects-upper.npy")
#     quad_centroids_np = np.load(folder_system/"output-centroids.npy")
#     n_frames = len(quad_centroids_np)

#     defect_sizes = []
#     defect_size_counts = []
#     defect_size_counts_log_transformed = []
#     bins_below50 = np.arange(0, 50, 2) #each bin is 2A
#     bins_above50 = np.arange(50, 451, 400)
#     bins = np.concatenate([bins_below50, bins_above50])
#     for frame in range(0, n_frames):
#         defects_counts_frame =  defect_counts[frame]

#         #all defect mask
#         masked_all_defects_num = np.where(
#             defects_counts_frame < 1,
#             1,
#             0
#         )

#         centroids_frame =  quad_centroids_np[frame, :, 0:1]

#         #restrain to region of interest 
#         centered_centroids = centroids_frame  - np.mean(centroids_frame)
#         centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()
#         restricted_region_grids_with_defects =  np.where((centroid_X_mask == 1) & (masked_all_defects_num == 1),1,0)

#         #quantifying the size of the defect 
#         label_all, nfeat_all = scipy_ndimage.label( restricted_region_grids_with_defects)
#         for i in range(1, nfeat_all + 1):
#             defect_size = np.count_nonzero(label_all == i)
#             defect_sizes.append(defect_size)
    
#     print(composition, "max defect size:",np.max(defect_sizes) )
#     #bin the data
#     data1, _ = np.histogram(
#         defect_sizes, bins = bins
#     )
#     bin_centers = 0.5 * (bins[1:] + bins[:-1]) 

#     #log transform
#     log_transform = np.log(data1)

#     bin_widths = np.diff(bins)

    
#     system_defect_sizes_log_transformed_count[composition] = log_transform


# with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
#     #fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.5, subheight=1.200)
#     fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=3)
#     alpha = 1
#     #colors_presets = [ "#5db7fc","#f8c362ff"]
#     colors_presets = [ "#d2826cff","#aaa9a9ff"]
#     i = 0
#     for composition in lipid_comp:
#         color = colors_presets[i]
#         ax.bar(bin_centers, system_defect_sizes_log_transformed_count[composition], alpha = alpha, width = bin_widths, label = composition, color = color)
#         ax.set_ylabel("log(defect-counts)")
#         ax.xaxis.set_major_locator(MultipleLocator(10))
#         ax.set_xlabel("Defect Area (Å²)")
#         alpha = alpha 
#         ax.set_xlim(0,55)
#         i += 1
#     plt.legend()
#     plt.savefig("logdefect-barplot-fig-peak-PC-ext.pdf", dpi = 300)
#     plt.savefig("logdefect-barplot-fig-peak-PC-ext.png", dpi = 300)
    










# # violin plot (overlaid)
# # log transform the counts 

# with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):

#     fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=3)


#     for composition, (areas, log_counts) in system_defect_sizes_log_transformed_count.items():
#         ax.scatter(areas, log_counts, label=composition, alpha=0.7)

#     ax.set_ylabel("log(Defect Count)")
#     ax.set_xlabel("size of defect")
#     ax.set_title("Distribution of Log-Transformed Defect Counts")
#     plt.legend()
#     plt.savefig(
#         f"new-violin-plot-defect-sizes-{configurations}-restriction{restriction}-PC.png",
#         dpi=300
#     )
#     plt.close()





# with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):

#     fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=3)

#     data = []
#     labels = []

#     for composition in lipid_compositions:
#         data.append(log_transformed_defect_size[composition])
#         labels.append(composition)

#     position = 1  # single x location for stacking
#     colors = ["#4C72B0", "#DDDB52"]
#     legend_handles = []

#     for i, dataset in enumerate(data):
#         vp = ax.violinplot(
#             dataset,
#             positions=[position],
#             showextrema=False,
#             widths=0.8
#         )

#         body = vp["bodies"][0]
#         body.set_facecolor(colors[i])
#         body.set_edgecolor("black")
#         body.set_alpha(1)
#         body.set_label(labels[i])   # ← key line

#     ax.set_xlabel("System")
#     ax.set_ylabel(r"Defect size ($\AA^2$)")
#     ax.legend(title="Lipid Composition", frameon=False)

#     plt.savefig(
#         f"new-violin-plot-defect-sizes-{configurations}-restriction{restriction}-PC.png",
#         dpi=300
#     )
#     plt.close()



# def total_num_defect_restricted_Z(data,shape, centroids_frame):

#     percent_cutoff = .2

#     reshaped_defects =  data.reshape(shape) #reshape the data so it follows the organization of the grids (x,y)
    

#     #create a mask for the peak based on Z height 
#     reshaped_Z_coord_centriods = centroids_frame[:,2] #reshape the z coord of the centroids to follow the organization of the grids
#     peak_Z = reshaped_Z_coord_centriods.max()
#     # Z_mask = (reshaped_Z_coord_centriods >= peak_Z - 50)
#     n_top = int(percent_cutoff * reshaped_Z_coord_centriods.size)
#     top_idx = np.argsort(peak_Z - reshaped_Z_coord_centriods, axis=None)[:n_top]
#     Z_mask = np.zeros_like(reshaped_Z_coord_centriods, dtype=bool)
#     Z_mask[top_idx] = True
#     total_grids= sum(Z_mask)
#     restricted_region_defects = reshaped_defects[Z_mask]


#     #compute average number of defects for region of interest
#     total_num_defects_restricted = np.sum(restricted_region_defects)

#     percent_defect_coverage = total_num_defects_restricted/total_grids

#     return total_num_defects_restricted,  total_grids,  percent_defect_coverage 