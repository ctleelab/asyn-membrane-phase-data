#outline of code
#1) create map 
#2) average map over frames
#3) compress the y dimension down 


#map defects in grid space, there will an issue as the same grid doesn't represent the same area in the membrane through each frame


#plot defects
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from scipy import ndimage as scipy_ndimage


configuration = "strain2"
defect_cut_off = 10
time = "100ps"
defect_type = ["deep", "shallow", "all"]


lipid_compositions = {
    "DOPC": [0, 0],
    "DPPC": [0, 1],
    "DOPC-DOPA": [1, 0],
    "DPPC-DPPA": [1, 1],
    "DOPC-DOPS": [2, 0],
    "DPPC-DPPS": [2, 1],
}


analysis_path = f"/home/casakurai/scratch/asyn-phase-binding-data/analysis/defect-data-{time}"
analysis_path = Path(analysis_path)
analysis_defect_path = analysis_path/f"defect-cut-off-{defect_cut_off}A"
compiled_data = Path(f"/home/casakurai/scratch/asyn-phase-binding-data/Figures/defect-data/{configuration}/defect-cut-off-{defect_cut_off}A")
compiled_data.mkdir(exist_ok=True)

#dictionary containing the systems for spatial weight plot 
systems_dict_shallow_spatial_weight_along_x = {} 
systems_dict_deep_spatial_weight_along_x = {}
systems_dict_all_spatial_weight_along_x = {}

systems_dict_avg_defect_shallow_per_grid = {} 
systems_dict_avg_defect_deep_per_grid = {}
systems_dict_avg_defect_all_per_grid  = {}

systems_dict_height = {}

for composition, position in lipid_compositions.items():
    folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
    defect_counts = np.load(folder_system/"output-defects-upper.npy")
    quad_centroids_np = np.load(folder_system/"output-centroids.npy")
    shape = np.load(folder_system/"output-shape.npy")
    X = np.load(folder_system/"output-Xshape.npy")
    Y = np.load(folder_system/"output-Yshape.npy")

    n_frames = len(quad_centroids_np)

    
    deep_defect_counts_per_frame_along_x = np.empty((n_frames, shape[0]))
    shallow_defect_counts_per_frame_along_x = np.empty((n_frames,shape[0]))
    all_defect_counts_per_frame_along_x = np.empty((n_frames,shape[0]))

    deep_defect_counts_per_frame_per_grid = np.empty((n_frames, shape[0], shape[1]))
    shallow_defect_counts_per_frame_per_grid = np.empty((n_frames,shape[0], shape[1]))
    all_defect_counts_per_frame_per_grid = np.empty((n_frames,shape[0], shape[1]))

    z_coord_per_frame = np.empty((n_frames,shape[0]))

    for frame in np.arange(0, n_frames):
        defect_counts_frame = defect_counts[frame]
        centroid_frame = quad_centroids_np[frame]

        z_coord_centroid_frame = centroid_frame[:,2]

        #create masks to denote presence of absence of defect

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

        #reshape mask to represent shape of grids 
        #gets a value per grid
        reshaped_shallow_defects =  masked_shallow_defects_num.reshape(shape)
        reshaped_deep_defects =  masked_deep_defects_num.reshape(shape)
        reshaped_all_defects =  masked_all_defects_num.reshape(shape)

        #sum down y, to get a count per grid in the x-dimension 
        sum_shallow_defects_along_y = np.sum(reshaped_shallow_defects, axis = 1)
        sum_deep_defects_along_y = np.sum(reshaped_deep_defects, axis = 1)
        sum_all_defects_along_y = np.sum(reshaped_all_defects, axis = 1)


        #save the sum per grid of defect counts along the x-axis per frame
        shallow_defect_counts_per_frame_along_x[frame, :] = sum_shallow_defects_along_y
        deep_defect_counts_per_frame_along_x[frame, :] = sum_deep_defects_along_y
        all_defect_counts_per_frame_along_x[frame, :] = sum_all_defects_along_y

        #save the grid value per frame
        shallow_defect_counts_per_frame_per_grid[frame, :] = reshaped_shallow_defects
        deep_defect_counts_per_frame_per_grid[frame, :] = reshaped_deep_defects
        all_defect_counts_per_frame_per_grid[frame, :] = reshaped_all_defects


        #height calculations
        reshaped_z_coord = z_coord_centroid_frame.reshape(shape)
        sum_z_coord_along_y = np.sum(reshaped_z_coord, axis = 1)
        z_coord_per_frame[frame, :] = sum_z_coord_along_y


    #sum of counts at a particular x grid across all frames 
    sum_deep_defect_across_frames_along_x =  np.sum(deep_defect_counts_per_frame_along_x, axis = 0)
    sum_shallow_defect_across_frames_along_x =  np.sum(shallow_defect_counts_per_frame_along_x, axis = 0)
    sum_all_defect_across_frames_along_x =  np.sum(all_defect_counts_per_frame_along_x, axis = 0)

    #divide by the number of frames 
    #average count at a particular grid 
    avg_deep_defect_across_frames_along_x =  sum_deep_defect_across_frames_along_x/n_frames
    avg_shallow_defect_across_frames_along_x =  sum_shallow_defect_across_frames_along_x/n_frames
    avg_all_defect_across_frames_along_x =  sum_all_defect_across_frames_along_x/n_frames


    #total average defects in one slice of y across all frames 
    sum_total_avg_deep_defects_across_frames_along_x = np.sum(avg_deep_defect_across_frames_along_x, axis = 0 )
    sum_total_avg_shallow_defects_across_frames_along_x = np.sum(avg_shallow_defect_across_frames_along_x, axis = 0)
    sum_total_avg_all_defects_across_frames_along_x = np.sum(avg_all_defect_across_frames_along_x, axis = 0)
    
    #percentage of defects in particular region of the membrane
    percent_deep_defects_x_dimension_along_x = np.divide(avg_deep_defect_across_frames_along_x, sum_total_avg_deep_defects_across_frames_along_x)
    percent_shallow_defects_x_dimension_along_x = np.divide(avg_shallow_defect_across_frames_along_x, sum_total_avg_shallow_defects_across_frames_along_x)
    percent_all_defects_x_dimension_along_x = np.divide(avg_all_defect_across_frames_along_x, sum_total_avg_all_defects_across_frames_along_x)

    #handeling to get average z-coord per x 
    avg_z_coord_across_frames = np.mean(z_coord_per_frame, axis = 0)

    #add to system dict
    systems_dict_shallow_spatial_weight_along_x[composition] = percent_shallow_defects_x_dimension_along_x
    systems_dict_deep_spatial_weight_along_x[composition] = percent_deep_defects_x_dimension_along_x
    systems_dict_all_spatial_weight_along_x[composition] = percent_all_defects_x_dimension_along_x

    #calculate the mean value of defect count add to system dict 
    systems_dict_avg_defect_shallow_per_grid[composition] = np.mean(shallow_defect_counts_per_frame_per_grid, axis = 0)
    systems_dict_avg_defect_deep_per_grid[composition] = np.mean(deep_defect_counts_per_frame_per_grid, axis = 0)
    systems_dict_avg_defect_all_per_grid[composition] = np.mean(all_defect_counts_per_frame_per_grid, axis = 0)

    systems_dict_height[composition] = avg_z_coord_across_frames


#percentage of defects in particular region of the membrane# 
#tells you where defects tend to be located along x
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
            ax.set_xlabel("X (Å)")
            ax.set_ylabel("Fraction of total defect density")
            ax.set_xlim(-len(compiled_defect_counts)/2, len(compiled_defect_counts)/2)
            ax.set_ylim(0, .02)
            ax.set_title(composition, fontsize=16)
    plt.savefig(f"{type}-compiled-defects-averaged-frames-spatial-weight-{configuration}")
    plt.close()


#plot height of membranes 
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
        
        height_membrane = systems_dict_height[composition]


        #denotes where the graph will lie
        i, j = position
        ax = fig.add_subplot(gs[i,j])
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


#plot for heatmap, where each grid has a value representing the probability of having a defect
for type in defect_type:
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig = plt.figure(figsize=(10,6))
    
        gs = gridspec.GridSpec(
            3, 3,
            width_ratios=[1, 1, 0.05],  
            wspace=0.3,
            hspace=0.1
        )

        ax_array = np.empty((3, 2), dtype=object)

        for composition, position in lipid_compositions.items():

            folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
            if type == "shallow":
                compiled_defect_counts = systems_dict_avg_defect_shallow_per_grid[composition]
            if type == "deep":
                compiled_defect_counts = systems_dict_avg_defect_deep_per_grid[composition]
            if type == "all":
                compiled_defect_counts = systems_dict_avg_defect_all_per_grid[composition]

            
            shape = np.load(folder_system/"output-shape.npy")
            X = np.load(folder_system/"output-Xshape.npy")
            Y = np.load(folder_system/"output-Yshape.npy")

            X_centered = X - X.mean()


            # get subplot position
            i, j = position
            ax = fig.add_subplot(gs[i,j])
            ax_array[i, j] = ax  # pick the correct Axes to plot the system

            # plot
            sc = ax.pcolor(X_centered, Y, compiled_defect_counts, cmap="viridis", shading="auto", vmin = 0, vmax = .5)

            # formatting
            ax.set_aspect("equal", adjustable="datalim")
            ax.xaxis.set_major_locator(MultipleLocator(25))
            ax.set_xlabel("X (A)")
            ax.set_ylabel("Y (A)")
            ax.set_ylim(Y.min(), Y.max())
            ax.set_aspect('equal', adjustable='box') 
            ax.set_title(composition, fontsize=16)
            cax = fig.add_subplot(gs[:, 2])
            cbar = fig.colorbar(sc, cax=cax)
            cbar.set_label(f"Probability of finding a {type} defect")
    plt.savefig(f"{type}-compiled-defects-averaged-frames-heatmap-{configuration}")
    plt.close()




# #plot for defect size# 
# for type in defect_type:
#     with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
#         fig = plt.figure(figsize=(10,6))
    
#         gs = gridspec.GridSpec(
#             3, 3,
#             width_ratios=[1, 1, 0.05],  
#             wspace=0.3,
#             hspace=1
#         )

#         ax_array = np.empty((3, 2), dtype=object)

#         for composition, position in lipid_compositions.items():
#             folder_system = analysis_defect_path / f"{composition}-{configuration}-production-stripped_centered"
#             defect_counts = np.load(folder_system/"output-defects-upper.npy")

#             all_defect_sizes = []

#             for frame in np.arange(0, len(defect_counts)):
#                 frame_defects = defect_counts[frame]
#                 #msk of all defects
#                 masked_all_defects_num = np.where(
#                     frame_defects < 1,
#                     1,
#                     0
#                 )


#                 # LABEL CONNECTED COMPONENTS IN BINARY DEFECT MAP
#                 label_all, nfeat_all = scipy_ndimage.label(masked_all_defects_num)
                
#                 for i in range(1, nfeat_all + 1):
#                     mask = np.where(label_all == i)
#                     all_defect_sizes.append(np.count_nonzero(mask))
            

#             # #divide by the number of frames to get an average
#             # all_defect_sizes_np = np.array(all_defect_sizes)
#             # average_all_defect_sizes_np = all_defect_sizes_np/600 #divide by the number of frames 

#             #denotes where the graph will lie
#             i, j = position
#             ax = fig.add_subplot(gs[i,j])
#             ax_array[i, j] = ax  
            
#             ax.hist(all_defect_sizes, bins=50)

#             ax.set_xlabel("Defect size (A^2)")
#             ax.set_ylabel("Count")
#             ax.set_title(composition, fontsize=16)
#     plt.savefig(f"{type}-compiled-defect-size-averaged-frames-{configuration}")
#     plt.close()

