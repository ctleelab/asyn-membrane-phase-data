#plotting defects

#import numpy array
##defect count
##grid centroids

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph

#plot for defect over top 10% grid height 

systems = ["DOPC-strain2", "DPPC-strain2", "DOPC-DOPS-strain2", "DPPC-DPPS-strain2", "DOPC-DOPA-strain2", "DPPC-DPPA-strain2"]
lipid_comp = ["DOPC", "DPPC", "DOPC-DOPS", "DPPC-DPPS", "DOPC-DOPA", "DPPC-DPPA"]
systems_dict ={} #dictionary containing the systems 

def average_defect_peak_10(defect_count, triangle_centroids, system):
    # Load the numpy arrays
    defect_count_np = np.load(defect_count)
    triangle_centroids_np = np.load(triangle_centroids)

    #create a numpy array (number frames, 2)
    #2 values is (defects per peak, defect total area)
    list_defects = np.empty((len(triangle_centroids_np),2))

    #loop over frames 
    for indx in np.arange(0, len(triangle_centroids_np)):
        total_num_grids = triangle_centroids_np.shape[1]
        Z = triangle_centroids_np[indx][:,2]
        defect_count_grid = defect_count_np[indx]

        n_top = int(0.20 * Z.size)
        z_max = np.max(Z)

        # rank-based mask
        top_idx = np.argsort(z_max - Z)[:n_top]
        mask = np.zeros_like(Z, dtype=bool)
        mask[top_idx] = True

        #total # of grids with height in the top 10% 
        total_grid_80 = sum(mask)



        #pull out the number of those that have defects <1 
        defect_grid_80 = sum((mask) & (defect_count_grid < 1))
        defect_grid_total= sum((defect_count_grid < 1))

        #percent defect 
        defect_percent_area_peak = (defect_grid_80/total_grid_80)*100
        defect_percent_area_total = (defect_grid_total/total_num_grids)*100

        list_defects[indx] = [defect_percent_area_peak,defect_percent_area_total]

    systems_dict[system] = list_defects

for system in systems:
    time ="10"
    file_path = Path(f"/home/casakurai/scratch/packmem/{system}-production-stripped-centered-{time}ns")
    defect_count = "output-defect-count-per-frame-upper.npy"
    triangle_centroids = "output-triangle-centroids-per-frame-upper.npy"
    average_defect_peak_10(file_path/defect_count, file_path/triangle_centroids,system)

# Save the results
with open("all-systems-defects.txt", "w") as f:
    for key, value in systems_dict.items():
        f.write(f"{key}: {value}\n")

aggregated_data = {}
for key in systems_dict:
    # Pull out peak defects and total defects
    peak_defects = systems_dict[key][:, 0]
    total_defects = systems_dict[key][:, 1]
    average_peak_defects_over_frame = np.mean(peak_defects)
    average_total_defects_over_frame = np.mean(total_defects)
    std_peak_defects_over_frame = np.std(peak_defects)
    std_total_defects_over_frame = np.std(total_defects)

    aggregated_data[key] = [average_peak_defects_over_frame, average_total_defects_over_frame, std_peak_defects_over_frame, std_total_defects_over_frame]

# Plot the aggregated data
def bar_graph_system(data):
    systems = list(data.keys())
    avg_peak_defects = [data[sys][0] for sys in systems]
    avg_total_defects = [data[sys][1] for sys in systems]
    std_peak_defects = [data[sys][2] for sys in systems]
    std_total_defects = [data[sys][3] for sys in systems]

    x = np.arange(len(systems))  
    width = 0.35 
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=1.5)
        bar1 = ax.bar(x - width/2, avg_peak_defects, width, yerr=std_peak_defects, label='Peak Defects', capsize=5)
        bar2 = ax.bar(x + width/2, avg_total_defects, width, yerr=std_total_defects, label='Total Defects', capsize=5)

        ax.set_xlabel('Systems')
        ax.set_ylabel('Defect Percentage')
        ax.set_title('Average Defect across all trajectories')
        ax.set_xticks(x)
        ax.set_xticklabels(lipid_comp, rotation=45, ha='right')
        ax.legend()

        # Add value annotations on top of bars
        for bar_group in [bar1, bar2]:
            for bar in bar_group:
                height = bar.get_height()
                ax.annotate(f'{height:.1f}',
                            xy=(bar.get_x() + bar.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom')

        fig.tight_layout()
        plt.savefig("bar-plt-all-systems.png")
        plt.show()

bar_graph_system(aggregated_data)



