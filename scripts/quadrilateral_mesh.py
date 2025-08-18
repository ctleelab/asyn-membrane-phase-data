
import pickle
import numpy as np
import MDAnalysis as mda
from membrane_spectral_analysis.base import MembraneSpectralAnalysis
from tqdm.contrib.concurrent import process_map
from util import analysis_path

#Cailyns outline for code
#goal create a mesh from the PO4 beads in each leaflets
#determine the min dimension of the system. Since the system is NVT there should be no change in the box size
#determine the size of the bins
#both in x and y direction
#divide the min dimension by the bin size to get the number of bins 
#MembraneSpectralAnalysis 

#size of bins (2nm)
step = 2

path = analysis_path / "curvature_selection" / "mesh"

def mesh(args):
    gro, selection, traj, path = args
    u = mda.Universe(gro, str(traj) , continous=True)

    # Average box dimensions over trajectory
    avg_dim_x = np.mean([ts.dimensions[0] for ts in u.trajectory])
    avg_dim_y = np.mean([ts.dimensions[1] for ts in u.trajectory])

    # Convert step to Å (MDAnalysis uses Å)
    step_angstrom = step * 10.0  # step in nm

    # Determine number of bins
    num_bins_x = int(avg_dim_x / step_angstrom)
    num_bins_y = int(avg_dim_y / step_angstrom)

    # Create bin edges
    x_edges = np.linspace(0, avg_dim_x, num_bins_x + 1)
    y_edges = np.linspace(0, avg_dim_y, num_bins_y + 1)

    # Bin centers for plotting
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    X, Y = np.meshgrid(x_centers, y_centers, indexing="ij")

    # Checking z-range of phosphates
    ag = u.atoms.select_atoms("name PO4")
    z_coords = ag.positions[:, 2]
    print("min z:", np.min(z_coords))
    print("max z:", np.max(z_coords))
    print("z spread:", np.ptp(z_coords))

    # Run MembraneSpectralAnalysis
    mc = MembraneSpectralAnalysis(
        u,
        select=selection,
        n_x_bins=num_bins_x,
        n_y_bins=num_bins_y,
        x_range=(0, avg_dim_x),
        y_range=(0, avg_dim_y),
        wrap=True,          # enable wrapping to avoid edge artifacts
        interpolate=True
    ).run(verbose=True)

    # Save result
    with open(path / "membrane_curvature_smooth.pickle", "wb") as handle:
        pickle.dump(mc, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return X, Y, mc  # optionally return mesh and result for plotting

jobs = []
systems = [2]
lipid = "DPPC"
time = 460000


for sys in systems:
    #gro file
    gro = analysis_path / "curvature_selection"/"NVT"/f"system{sys}-8x8x25-{lipid}-{time}ps-NVT"/"equil"/"equilibration7.6.gro"
    #traj file
    traj = analysis_path / "curvature_selection"/"NVT"/f"system{sys}-8x8x25-{lipid}-{time}ps-NVT"/"equil"/"equilibration7.6.xtc"
    #selection
    selection = f"name PO4"
    inputs = (gro, selection, traj, path)
    jobs.append(inputs)

process_map(mesh, jobs, max_workers=64)




