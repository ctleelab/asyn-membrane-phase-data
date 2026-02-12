

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
from scipy import ndimage as scipy_ndimage


configurations = ["flat"]
defect_cut_off = 10
time = "100ps"
defect_type = ["deep", "shallow", "all"]
lower_limit_size_defect = 15
upper_limit_size_defect = 60

# lipid_compositions = {
#     "DOPC": [0, 0],
#     "DPPC": [0, 1],
#     "DOPC-DOPA": [1, 0],
#     "DPPC-DPPA": [1, 1],
#     "DOPC-DOPS": [2, 0],
#     "DPPC-DPPS": [2, 1]
# }

lipid_compositions = {
    "DOPC-DOPA": [1, 0],
    "DPPC-DPPA": [1, 1]
}



lipid_comp = list(lipid_compositions.keys())

analysis_path = f"/home/casakurai/scratch/asyn-phase-binding-data/analysis/defect-data-{time}"
# analysis_path = Path(analysis_path)
analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")

# # #############################################################################
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




# # def bar_graph_system_peak_grouped(data,restriction,configurations):
# #     systems = list(data.keys())
# #     n_configs = len(configurations)

# #     #set labels for system
# #     x = np.arange(len(systems))  

# #     width = 0.35 
# #     system_spacing = 4   # separates systems
# #     config_spacing = 1  # keeps configs tight
    
# #     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
# #         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.5, subheight=1.25)
        

# #         #plot each configuration
# #         for i, config in enumerate(configurations):
# #             avg_values = [data[sys][i][0] for sys in systems]
# #             std_values = [data[sys][i][1] for sys in systems]
# #             group_spacing = 1.2 
# #             offsets = (i - (n_configs - 1) / 2) * width * config_spacing

# #             ax.bar(x + offsets,  avg_values, width ,yerr = std_values, label=config, capsize=5)
# #         ax.legend()
# #         ax.set_ylabel('Surface Covered by Defects (%)')
# #         ax.set_xticks(x)
# #         #ax.set_ylim(0,10)

# #         ax.set_xticklabels(lipid_comp, rotation=45, ha='right')

# #         plt.savefig(f"defect-coverage-bar-plt-avg-frames-allconfigurations-restriction{restriction}-PC.pdf")

# # # # def total_num_defect_restricted_Z(data,shape, centroids_frame):

# # # #     percent_cutoff = .2

# # # #     reshaped_defects =  data.reshape(shape) #reshape the data so it follows the organization of the grids (x,y)
    

# # # #     #create a mask for the peak based on Z height 
# # # #     reshaped_Z_coord_centriods = centroids_frame[:,2] #reshape the z coord of the centroids to follow the organization of the grids
# # # #     peak_Z = reshaped_Z_coord_centriods.max()
# # # #     # Z_mask = (reshaped_Z_coord_centriods >= peak_Z - 50)
# # # #     n_top = int(percent_cutoff * reshaped_Z_coord_centriods.size)
# # # #     top_idx = np.argsort(peak_Z - reshaped_Z_coord_centriods, axis=None)[:n_top]
# # # #     Z_mask = np.zeros_like(reshaped_Z_coord_centriods, dtype=bool)
# # # #     Z_mask[top_idx] = True
# # # #     total_grids= sum(Z_mask)
# # # #     restricted_region_defects = reshaped_defects[Z_mask]


# # # #     #compute average number of defects for region of interest
# # # #     total_num_defects_restricted = np.sum(restricted_region_defects)

# # # #     percent_defect_coverage = total_num_defects_restricted/total_grids

# # # #     return total_num_defects_restricted,  total_grids,  percent_defect_coverage 




def total_num_defect_restricted_X(defect_data,shape,centroids,restriction,configuration):
    
    #center of grids at 0 
    centered_centroids = centroids - np.mean(centroids)

    centroid_X_mask = (np.where((centered_centroids >= -restriction) & (centered_centroids <= restriction),1,0)).squeeze()

    restricted_region_grids_with_defects =  np.where((centroid_X_mask == 1) & (defect_data == 1),1,0)

    total_grids = len(defect_data)


    defect_grids = np.sum(restricted_region_grids_with_defects)



    percent_defect_coverage = (defect_grids/total_grids)*100

    # print("percent coverage:", percent_defect_coverage)

    #return total_num_defects_restricted_X, total_grids_restricted_X
    return defect_grids,  total_grids,  percent_defect_coverage 


# # ####Surface covered by defects ####
# # restrictions = [25]
# # defect_restricted_avg_system = {} 
# # for composition in lipid_comp:
# #     configurations_types = {}
# #     for configuration in configurations: 
# #         folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
# #         defect_counts = np.load(folder_system/"output-defects-upper.npy")
# #         quad_centroids_np = np.load(folder_system/"output-centroids.npy")
# #         shape = np.load(folder_system/"output-shape.npy")
# #         X = np.load(folder_system/"output-Xshape.npy")
# #         Y = np.load(folder_system/"output-Yshape.npy")

# #         n_frames = len(quad_centroids_np)

# #         # perecent_cut_defects = np.empty(6, n_frames)

# #         # for percent_cutoff in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:

# #         perecent_cut_defects = {}
# #         for restriction in restrictions: 

# #             total_defects_restricted_region_per_frame = np.empty(n_frames)
# #             total_grids_restricted_region_per_frame = np.empty(n_frames)
# #             percentage_defect_coverage_per_frame = np.empty(n_frames)

# #             for frame in np.arange(0, n_frames):


# #                 defect_counts_frame = defect_counts[frame]
# #                 centroids_frame =  quad_centroids_np[frame, :, 0:1]

# #                 #all defect 
# #                 masked_all_defects_num = np.where(
# #                     defect_counts_frame < 1,
# #                     1,
# #                     0
# #                 )

# #                 total_per_frame, total_grids, percent_defect_coverage = total_num_defect_restricted_X(masked_all_defects_num,shape,centroids_frame,restriction,configuration)

# #                 total_defects_restricted_region_per_frame[frame] = total_per_frame

# #                 total_grids_restricted_region_per_frame [frame] = total_grids

# #                 percentage_defect_coverage_per_frame[frame] = percent_defect_coverage

# #             average_defect_coverage = np.mean(percentage_defect_coverage_per_frame)
# #             std_defect_coverage = np.std(percentage_defect_coverage_per_frame)

# #             perecent_cut_defects[restriction]  = [average_defect_coverage,std_defect_coverage]
# #         configurations_types[configuration] = perecent_cut_defects
# #     defect_restricted_avg_system[composition]= configurations_types
    


# # restriction = 25
# # data_grouped = {}

# # for system in lipid_comp:
# #     data_grouped[system] = []
# #     for config in configurations:
# #         avg, std = defect_restricted_avg_system[system][config][restriction]
# #         data_grouped[system].append([avg, std])

# #     #bar_graph_system_peak(temp_dict, restriction)

# # bar_graph_system_peak_grouped(data_grouped, restriction, configurations)


# # ###########################################################
# # ###calculate the size of the defects in a defined region### 
# # ###########################################################

# restriction = 200
# system_defect_sizes = {}
# for composition, position in lipid_compositions.items():
#     folder_system = analysis_defect_path / f"{composition}-{configurations}-production-stripped_centered"
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


#         label_all, nfeat_all = scipy_ndimage.label( restricted_region_grids_with_defects)

#         for i in range(1, nfeat_all + 1):
#             defect_size = np.count_nonzero(label_all == i)
#             all_defect_sizes.append(defect_size)

#     system_defect_sizes[composition] = all_defect_sizes 

# # ############################################
# # ###calculate the maximum size of a defect### 
# # ############################################

# largest_defect = 0
# restriction = 200
# system_largest_defect_sizes = {}
# for composition, position in lipid_compositions.items():
#     largest_defect = 0
#     folder_system = analysis_defect_path / f"{composition}-{configurations}-production-stripped_centered"
#     defect_counts = np.load(folder_system/"output-defects-upper.npy")
#     quad_centroids_np = np.load(folder_system/"output-centroids.npy")

#     all_defect_sizes = [] #a list that states each defect size
#     n_frames = len(quad_centroids_np)


#     for frame in np.arange(0, n_frames):
#         defects_counts_frame =  defect_counts[frame]

#         #all defect 
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


#         label_all, nfeat_all = scipy_ndimage.label( restricted_region_grids_with_defects)

        
#         for i in range(1, nfeat_all + 1):
#             defect_size = np.count_nonzero(label_all == i)
#             if defect_size > largest_defect:
#                 largest_defect = defect_size
                

#     system_largest_defect_sizes[composition] = largest_defect 
# print(system_largest_defect_sizes)


# ##########################################
# ###bar plot avg defect size over frames### 
# ##########################################
# def bar_graph_system(data):
#     systems = list(data.keys())
#     list_of_defect_sizes  = [data[sys] for sys in systems]
#     avg_defects_size = [np.mean(sizes) for sizes in list_of_defect_sizes]
#     std_defect_size = [np.std(sizes) for sizes in list_of_defect_sizes]



#     #set labels for system
#     x = np.arange(len(systems))  

#     width = 0.35 
    
#     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
#         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=3)
#         #bar1 = ax.bar(x - width/2, avg_peak_defects, width, yerr=std_peak_defects, label='Peak Defects', capsize=5)
#         bar2 = ax.bar(x + width/2,  avg_defects_size, width, yerr = std_defect_size, label='Total Defects', capsize=5)

#         ax.set_xlabel('Systems')
#         ax.set_ylabel('Average defect size A^2')
#         ax.set_xticks(x)

#         ax.set_xticklabels(lipid_comp, rotation=45, ha='right')
#         ax.legend()

#         # Add value annotations on top of bars
#         for bar_group in [bar2]:
#             for bar in bar_group:
#                 height = bar.get_height()
#                 ax.annotate(f'{height:.1f}',
#                             xy=(bar.get_x() + bar.get_width() / 2, height),
#                             xytext=(0, 3), 
#                             textcoords="offset points",
#                             ha='center', va='bottom')


#         fig.tight_layout()
#         plt.savefig(f"defect-size-bar-plt-avg-frames-{configurations}.png")
#         plt.show()

# # bar_graph_system(system_defect_sizes )



# # ###################################
# # ###distrubution of defect sizes ### 
# # ###################################
# # #area under curve is 1, so shape of defect histograms are comparable

# system_defect_sizes = {}
# with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
#     fig,ax = ph.fixed_size_subplots(1, 1, subwidth=1, subheight=1)
#     for composition, position in lipid_compositions.items():
#         folder_system = analysis_defect_path / f"{composition}-{configurations}-production-stripped_centered"
#         defect_counts = np.load(folder_system/"output-defects-upper.npy")

#         all_defect_sizes = [] #a list that states has each defect size 
#         bins = np.arange(0, 400, 1) #define bins so that they are equal between systems

#         for frame_defects in defect_counts:
#             masked_all_defects_num = np.where((frame_defects < 1), 1, 0)

#             label_all, nfeat_all = scipy_ndimage.label(masked_all_defects_num)

#             for i in range(1, nfeat_all + 1):
#                 defect_size = np.count_nonzero(label_all == i)
#                 all_defect_sizes.append(defect_size)
#         system_defect_sizes[composition] = all_defect_sizes 


#         pdf, bins = np.histogram(all_defect_sizes, bins=bins, density=True)
#         bin_centers = 0.5 * (bins[1:] + bins[:-1])


#         ax.bar(bin_centers, pdf, label = composition)

# ax.set_ylabel("?")
# ax.set_xlim(0,55)
# ax.set_xlabel("Defect Area (A^2)")
# ax.legend(loc = "upper right", fontsize = "small")
# plt.tight_layout()
# plt.savefig(f"all-probability-{composition}-defect-size-{configurations}")
# plt.close()


# ########################################
# ###convergence of defect % coverage  ### 
# ########################################

restrictions = [25]
defect_restricted_system = {} 
defect_restricted_system_per_frame = {} 
for composition in lipid_comp:
    configurations_types = {}
    avg_defect_per_frame ={}
    for configuration in configurations: 
        folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
        defect_counts = np.load(folder_system/"output-defects-upper.npy")
        quad_centroids_np = np.load(folder_system/"output-centroids.npy")
        shape = np.load(folder_system/"output-shape.npy")
        X = np.load(folder_system/"output-Xshape.npy")
        Y = np.load(folder_system/"output-Yshape.npy")

        n_frames = len(quad_centroids_np)

        # perecent_cut_defects = np.empty(6, n_frames)

        # for percent_cutoff in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:

        perecent_cut_defects = {}
        perecent_cut_defects_per_frame ={}
        for restriction in restrictions: 

            total_defects_restricted_region_per_frame = np.empty(n_frames)
            total_grids_restricted_region_per_frame = np.empty(n_frames)
            percentage_defect_coverage_per_frame = np.empty(n_frames)
        
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

            perecent_cut_defects[restriction]  = [average_defect_coverage,std_defect_coverage]
            perecent_cut_defects_per_frame[restriction] = [percentage_defect_coverage_per_frame]
            
        configurations_types[configuration] = perecent_cut_defects
        avg_defect_per_frame[configuration] = perecent_cut_defects_per_frame
    defect_restricted_system[composition ]= configurations_types
    defect_restricted_system_per_frame[composition] = avg_defect_per_frame

#scatter plot for defect_restricted_system_per_frame
plt.figure(figsize=(10, 6))

for composition in defect_restricted_system_per_frame:
    for configuration in defect_restricted_system_per_frame[composition]:
        # Extract the per-frame defect coverage array
        per_frame_array = defect_restricted_system_per_frame[composition][configuration][restriction][0]  # [0] because you stored as list

        n_frames = len(per_frame_array)
        frames = np.arange(n_frames)

        # Scatter plot
        plt.scatter(frames, per_frame_array, label=f"{composition}-{configuration}", alpha=0.6)

plt.xlabel("Frame")
plt.ylabel("Percent defect coverage")
plt.title(f"Defect coverage per frame (restriction={restriction})")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig(f"frame-avg-defect-coverage{configurations}.png")