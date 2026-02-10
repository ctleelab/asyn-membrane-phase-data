
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, NullFormatter
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from scipy import ndimage as scipy_ndimage
from scipy.stats import linregress

from tqdm.auto import tqdm

#this code relies on the grid all being uniform, which is an assumption from the PackMem-leelab

configuration = "strain2" 
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
    "DOPC": [0, 0],
    "DPPC": [0, 1]}


lipid_comp = list(lipid_compositions.keys())

analysis_path = f"/home/casakurai/scratch/asyn-phase-binding-data/analysis/defect-data-{time}"
# analysis_path = Path(analysis_path)
analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
compiled_data = Path(f"/home/casakurai/scratch/asyn-phase-binding-data/Figures/defect-data-{time}/{configuration}/defect-cut-off-{defect_cut_off}A")
compiled_data.mkdir(exist_ok=True)

#dictionary containing the systems for spatial weight plot 
systems_dict_shallow_spatial_weight_along_x = {} 
systems_dict_deep_spatial_weight_along_x = {}
systems_dict_all_spatial_weight_along_x = {}

systems_dict_avg_defect_shallow_per_grid = {} 
systems_dict_avg_defect_deep_per_grid = {}
systems_dict_avg_defect_all_per_grid  = {}

systems_dict_height = {}

systemY_cutout_Y_artifacts ={}
systemX_cutout_Y_artifacts ={}

for composition, position in lipid_compositions.items():
    folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
    defect_counts = np.load(folder_system/"output-defects-upper.npy")
    quad_centroids_np = np.load(folder_system/"output-centroids.npy")
    shape = np.load(folder_system/"output-shape.npy")
    X = np.load(folder_system/"output-Xshape.npy")
    Y = np.load(folder_system/"output-Yshape.npy")

    n_frames = len(quad_centroids_np)
    cutsize = 8

    deep_defect_counts_per_frame_along_x = np.empty((n_frames, shape[0]))
    shallow_defect_counts_per_frame_along_x = np.empty((n_frames,shape[0]))
    all_defect_counts_per_frame_along_x = np.empty((n_frames,shape[0]))

    deep_defect_counts_per_frame_per_grid = np.empty((n_frames, shape[0], shape[1]-cutsize*2))
    shallow_defect_counts_per_frame_per_grid = np.empty((n_frames,shape[0], shape[1]-cutsize*2))
    all_defect_counts_per_frame_per_grid = np.empty((n_frames,shape[0], shape[1]-cutsize*2))

    z_coord_per_frame_along_x = np.empty((n_frames,shape[0]))

    for frame in tqdm(np.arange(0, n_frames)):
        defect_counts_frame = defect_counts[frame]
        centroid_frame = quad_centroids_np[frame]

        z_coord_centroid_frame = centroid_frame[:,2]

        # Add a value of 1 to the grid where the Boolean is true
        masked_shallow_defects_num = np.where(
            (defect_counts_frame >  0) & (defect_counts_frame < 1),
            1,
            0,
        )

        # Add a value of 1 to the grid where the Boolean is true
        masked_deep_defects_num = np.where(
            (defect_counts_frame ==  0),
            1,
            0,
        )

        #all defect 
        masked_all_defects_num = np.where(
            defect_counts_frame < 1,
            1,
            0
        )
    
        def cut_out_y_artifacts(arr, begin, end):
            out = np.hstack((arr[:, :begin], arr[:, end:,]))
            deleted = arr[:, begin:end]    
            return out, deleted 

        #reshape mask to represent shape of grids 
        #gets a value per grid
        reshaped_shallow_defects =  masked_shallow_defects_num.reshape(shape)
        reshaped_deep_defects =  masked_deep_defects_num.reshape(shape)
        reshaped_all_defects =  masked_all_defects_num.reshape(shape)

        center = X.shape[1] // 2

        min_index = center - cutsize
        max_index = center + cutsize
        # print(shape, center, min_index, max_index)

        reshaped_shallow_defects, _ = cut_out_y_artifacts(reshaped_shallow_defects, min_index, max_index)
        reshaped_deep_defects, _ = cut_out_y_artifacts(reshaped_deep_defects, min_index, max_index)
        reshaped_all_defects, _ = cut_out_y_artifacts(reshaped_all_defects, min_index, max_index)


        #new shape of array 
        new_shape = reshaped_shallow_defects.shape
        
        # print(stripped_reshaped_shallow_defects.shape, reshaped_shallow_defects.shape)
        # print(np.isclose(stripped_reshaped_shallow_defects[min_index-1], reshaped_shallow_defects[min_index-1]))
        # print(np.isclose(stripped_reshaped_shallow_defects[min_index], reshaped_shallow_defects[max_index]))
    

        #sum down y, to get a count per grid in the x-dimension 
        shallow_defect_counts_per_frame_along_x[frame] = np.sum(reshaped_shallow_defects, axis = 1)
        deep_defect_counts_per_frame_along_x[frame] = np.sum(reshaped_deep_defects, axis = 1)
        all_defect_counts_per_frame_along_x[frame] = np.sum(reshaped_all_defects, axis = 1)

        #save the grid value per frame
        shallow_defect_counts_per_frame_per_grid[frame] = reshaped_shallow_defects
        deep_defect_counts_per_frame_per_grid[frame] = reshaped_deep_defects
        all_defect_counts_per_frame_per_grid[frame] = reshaped_all_defects

        #height calculations
        reshaped_z_coord = z_coord_centroid_frame.reshape(shape)
        z_coord_per_frame_along_x[frame] = np.sum(reshaped_z_coord, axis = 1)


    #sum of counts at a particular x grid across all frames 
    sum_deep_defect_across_frames_along_x =  np.sum(deep_defect_counts_per_frame_along_x, axis = 0)
    sum_shallow_defect_across_frames_along_x =  np.sum(shallow_defect_counts_per_frame_along_x, axis = 0)
    sum_all_defect_across_frames_along_x =  np.sum(all_defect_counts_per_frame_along_x, axis = 0)

    sum_height_across_frames_along_x = np.sum(z_coord_per_frame_along_x, axis = 0)


    # #divide by the total number of grids (number of grids in y-dim * number of frames)
    total_grids = new_shape[1] * n_frames


    # #percentage of defects in particular region of the membrane
    systems_dict_deep_spatial_weight_along_x[composition]= np.divide(sum_deep_defect_across_frames_along_x, total_grids)
    systems_dict_shallow_spatial_weight_along_x[composition]  = np.divide(sum_shallow_defect_across_frames_along_x, total_grids)
    systems_dict_all_spatial_weight_along_x[composition] = np.divide(sum_all_defect_across_frames_along_x, total_grids)

    #handeling to get average z-coord per x 
    avg_z_coord_across_frames = np.divide(sum_height_across_frames_along_x, total_grids)

    #calculate the mean value of defect count add to system dict 
    systems_dict_avg_defect_shallow_per_grid[composition] = np.mean(shallow_defect_counts_per_frame_per_grid, axis = 0)
    systems_dict_avg_defect_deep_per_grid[composition] = np.mean(deep_defect_counts_per_frame_per_grid, axis = 0)
    systems_dict_avg_defect_all_per_grid[composition] = np.mean(all_defect_counts_per_frame_per_grid, axis = 0)

    min_avg_z_coord_across_frames = np.min(avg_z_coord_across_frames)
    max_avg_z_coord_across_frames = np.max(avg_z_coord_across_frames)

    avg_z_coord_across_frames_normalized = (avg_z_coord_across_frames - min_avg_z_coord_across_frames)/10 #convert to nm

    systems_dict_height[composition] = avg_z_coord_across_frames_normalized 
 
    kept_Y, deleted_Y  = cut_out_y_artifacts(Y,min_index, max_index)
    kept_X, deleted_X = cut_out_y_artifacts(X,min_index, max_index)

    #need to adjust the values in the Y and X to account for the deleted rows of data
    adjusted_kept_Y  = []
    for i in range (0, len(kept_Y[1])):
        number =  int(kept_Y[1][i])
        if number  > deleted_Y[0][0]:
            temp = int(number - cutsize*2)
            adjusted_kept_Y.append(temp)
        else:
            adjusted_kept_Y.append(number)
    
    
    #need to create the mesh with the cutout y-values
    x_axis = kept_Y.shape[0]
    X, Y = np.meshgrid(np.arange(x_axis), adjusted_kept_Y, indexing="ij")
    systemY_cutout_Y_artifacts[composition] = Y
    systemX_cutout_Y_artifacts[composition] = X



##############################################################################################
#########percentage of defects in particular region of the membrane in the x-dimension######### 
##############################################################################################

for type in defect_type:
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig = plt.figure(figsize=(10,6))

        gs = gridspec.GridSpec(
            3, 3,
            width_ratios=[1, 1, 0.05],  
            wspace=0.3,
            hspace=1
        )

        ax_array = np.empty((3, 2), dtype=object)

        for composition, position in lipid_compositions.items():

            #load the correct data 
            folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"

            if type == "shallow":
                compiled_defect_counts = systems_dict_shallow_spatial_weight_along_x [composition]
            if type == "deep":
                compiled_defect_counts = systems_dict_deep_spatial_weight_along_x[composition]
            if type == "all":
                compiled_defect_counts = systems_dict_all_spatial_weight_along_x [composition]


            #denotes where the graph will lie
            i, j = position
            ax = fig.add_subplot(gs[i,j])
            ax_array[i, j] = ax  
            
            
            x_idx = np.arange(len(compiled_defect_counts))

            x_idx_centered = x_idx - (len(compiled_defect_counts)/2)


            sc = ax.bar(x_idx_centered, 
                compiled_defect_counts,
                width=1,
                align="center"
            )

            # formatting
            ax.xaxis.set_major_locator(MultipleLocator(25))
            ax.set_ylim(0, .4)
            ax.set_xlabel("X (Å)")
            ax.set_ylabel("Fraction of total defect density")
            ax.set_xlim(-len(compiled_defect_counts)/2, len(compiled_defect_counts)/2)
            ax.set_title(composition, fontsize=16)
    plt.savefig(f"{type}-compiled-defects-averaged-frames-spatial-weight-{configuration}")
    plt.close()

########################################
#########height of the membrane######### 
########################################

with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
    fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=3)

    ax_array = np.empty((3, 2), dtype=object)

    for composition, position in lipid_compositions.items():

        #load the correct data 
        folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
        
        height_membrane = systems_dict_height[composition]


        #denotes where the graph will lie
        i, j = position
        ax_array[i, j] = ax  
        
        
        x_idx = np.arange(len(height_membrane))

        x_idx_centered = x_idx - (len(height_membrane)/2)

        sc = ax.bar(x_idx_centered, 
            height_membrane,
            width=1,
            align="center"
        )

        # formatting
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Z (Å)")
        ax.set_xlim(-len(height_membrane)/2, len(height_membrane)/2)
        ax.set_ylim(4000, 11000)
        ax.set_title(composition, fontsize=16)
plt.savefig(f"averaged-frames-height-{configuration}")
plt.close()

########################################################################################################################################
########combined graph of height of membrane and percentage of defects in particular region of the membrane in the x-dimension ######### 
########################################################################################################################################
for defect_type_name in defect_type:
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax = ph.fixed_size_subplots(1, 2, wmargin = .6, hmargin = 1,  subwidth=.95, subheight=1.1, rmargin_scale = 8)
        ax = np.array(ax).flatten()  # flatten in case it's 2D
        legend_handles = None
        legend_labels = None

        for idx, (composition, position) in enumerate(lipid_compositions.items()):
            # Load data
            height_membrane = systems_dict_height[composition]


            if defect_type_name == "shallow":
                compiled_defect_counts = systems_dict_shallow_spatial_weight_along_x[composition]
            elif defect_type_name == "deep":
                compiled_defect_counts = systems_dict_deep_spatial_weight_along_x[composition]
            else:
                compiled_defect_counts = systems_dict_all_spatial_weight_along_x[composition]

            x_idx = np.arange(len(height_membrane))
            x_idx_centered = (x_idx - (len(height_membrane)/2)) / 10   # convert to nm

            ax[idx].scatter(x_idx_centered, height_membrane, s=7, color="blue", label="Height")
            ax[idx].set_xlabel("X (nm)")
            ax[idx].set_ylabel("Height (nm)")
            ax[idx].set_xlim(-len(height_membrane)/20, len(height_membrane)/20)  # adjust for nm
            ax[idx].set_ylim(0, 15)  # nm
            ax[idx].set_title(f"{composition}", fontsize=12)  # nm
            ax[idx].xaxis.set_minor_locator(MultipleLocator(2.5))
            ax[idx].xaxis.set_major_locator(MultipleLocator(5))

            # Secondary y-axis (defect counts)

            ax2 = ax[idx].twinx()
            ax2.bar(x_idx_centered, compiled_defect_counts, width=.1, color="orange", alpha=0.5, label="Defect probability")
            ax2.set_ylabel("Defect probability")
            ax2.set_ylim(0,.5)
            

            if legend_handles is None:
                h1, l1 = ax[idx].get_legend_handles_labels()
                h2, l2 = ax2.get_legend_handles_labels()
                legend_handles = h1 + h2
                legend_labels = l1 + l2
        plt.legend(
            legend_handles,
            legend_labels,
            loc="center left",
            bbox_to_anchor=(1.2, 0.5),
            frameon=False,
            fontsize = 10
        )
        plt.subplots_adjust(right=0.82)
        
        plt.savefig(f"combined-defect-probability-and-height-{configuration}-{defect_type_name}-PC.png")
        plt.savefig(f"combined-defect-probability-and-height-{configuration}-{defect_type_name}-PC.pdf")
        plt.close()


#######################################################################################################
#########heatmap, where each grid has a value representing the probability of having a defect ######### 
#######################################################################################################
for type in defect_type:
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax = ph.fixed_size_subplots(1, 2, wmargin = .5, hmargin = 1,  subwidth=1.1, subheight=2)
        ax = np.array(ax).flatten() 

        for idx, (composition, position) in enumerate(lipid_compositions.items()):

            Y = systemY_cutout_Y_artifacts[composition]
            X = systemX_cutout_Y_artifacts[composition]
            

            folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
            if type == "shallow":
                compiled_defect_counts = systems_dict_avg_defect_shallow_per_grid[composition]
            if type == "deep":
                compiled_defect_counts = systems_dict_avg_defect_deep_per_grid[composition]
            if type == "all":
                compiled_defect_counts = systems_dict_avg_defect_all_per_grid[composition]
            
            print("compiled",compiled_defect_counts.shape)

            #need to reshape Y so it fits the data
            X_centered = X - X.mean()
            X_centered_nm = X_centered/10
            Y_nm = Y/10

            

            # plot
            ax[idx].pcolor(X_centered_nm ,Y_nm, compiled_defect_counts, cmap="viridis", shading="auto", edgecolor = "face", vmin = 0, vmax = .5)

            # formatting
            ax[idx].xaxis.set_major_locator(MultipleLocator(5))
            ax[idx].xaxis.set_minor_locator(MultipleLocator(2.5))
            ax[idx].yaxis.set_major_locator(MultipleLocator(2.5))
            ax[idx].set_xlabel("X (nm)")
            ax[idx].set_ylabel("Y (nm)")
            ax[idx].set_aspect('equal', adjustable='box') 
            ax[idx].set_title(composition, fontsize=16)
    plt.savefig(f"{type}-compiled-defects-averaged-frames-heatmap-{configuration}.pdf")
    plt.close()


# #####################################################################
# ###plot for the probability of finding a defect of a certain size### 
# #####################################################################

# #custom axis to match PackMem data
# major_probs = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]) 

# system_defect_sizes = {}
# with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
#     fig,ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=3)
#     for composition, position in lipid_compositions.items():
#         folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
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

#         #creating a histogram to represent the probability
#         counts, _ = np.histogram(
#             all_defect_sizes, bins = bins
#         )
#         probabilities = counts / counts.sum()
#         bin_centers = 0.5 * (bins[1:] + bins[:-1])

#         mask = probabilities > 0
#         areas_nonzero = bin_centers[mask]
#         log_probabilities = np.log(probabilities[mask])
        
#         #apply limits of defect size for fitting
#         masks_limits = (areas_nonzero >= lower_limit_size_defect)
#         areas_nonzero = areas_nonzero[masks_limits]
#         log_probabilities = log_probabilities[masks_limits]

#         #scatterplot 
#         ax.scatter(areas_nonzero, log_probabilities, s = 5, alpha =.7, label = composition)

#         #linear regression fit 
#         linear_regression = linregress(areas_nonzero, log_probabilities)

#         slope = linear_regression.slope
#         intercept = linear_regression.intercept
#         r2_value = linear_regression.rvalue**2 

#         xfit = np.linspace(min(all_defect_sizes), 100, 100)
#         yfit = slope * xfit + intercept 

# ax.plot(xfit, yfit, label = f"y={slope:.2f}x+{intercept:.2f}, R^2={r2_value:.2f}")
# ax.set_yticks(np.log(major_probs))
# ax.set_yticklabels([r"$10^{-1}$", r"$10^{-2}$", r"$10^{-3}$", r"$10^{-4}$", r"$10^{-5}$", r"$10^{-6}$"])
# ax.set_ylabel("ln(Probability)")
# ax.set_xlabel("Defect Area (A^2)")
# ax.set_ylim(np.log(1e-6), 0) 
# ax.legend(loc = "upper right", fontsize = "small")
# plt.tight_layout()
# plt.savefig(f"all-probability-{composition}-defect-size-{configuration}")
# plt.close()


# ######################################################
# ###violin scatterplot to visualize all defect sizes### 
# ######################################################

# def violin_scatter_plot(data):
#     systems = list(data.keys())
#     list_of_defect_sizes = [data[sys] for sys in systems]

#     x = np.arange(1, len(systems) + 1)

#     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
#         fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=1.5)

#         for i, sizes in enumerate(list_of_defect_sizes):
#             # jitter around x position
#             jitter = np.random.normal(loc=0, scale=0.05, size=len(sizes))
#             ax.scatter(
#                 np.full(len(sizes), x[i]) + jitter,
#                 sizes,
#                 s=1,
#                 alpha=0.6,
#             )

#         ax.set_xlabel("Systems")
#         ax.set_ylabel("Defect size (Å$^2$)")
#         ax.set_xticks(x)
#         #ax.set_ylim(15,100)
#         ax.set_xticklabels(lipid_comp, rotation=45, ha="right")

#         plt.savefig(f"defect-size-scatter-{configuration}.png", dpi=300)
#         plt.show()

# violin_scatter_plot(system_defect_sizes)



