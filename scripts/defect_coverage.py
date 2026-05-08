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


def total_num_defect_restricted_X_cutoff(
    defect_data, shape, centroids, restriction, defect_size_cutoff
):
    """
    Computed per frame
    """

    # reshapes to the original spatial (x,y)
    defect_grid = defect_data.reshape(shape)

    # center of grids at 0
    centered_centroids = centroids - np.mean(centroids)

    # creates a mask where only considered defects within the restricted area
    centroid_X_mask = (
        np.where(
            (centered_centroids >= -restriction) & (centered_centroids <= restriction),
            1,
            0,
        )
    ).squeeze()

    # reshapes the centroid_X_mask
    centroid_mask_grid = centroid_X_mask.reshape(shape)

    # combines the two masks (centroid_X_mask and defect_grid), so only get the data for defects within a certain range
    restricted_region_grids_with_defects = np.where(
        (centroid_mask_grid == 1) & (defect_grid == 1), 1, 0
    )

    # here handle the defect based on size
    all_defect_sizes_cutoff = []
    label_all, nfeat_all = scipy_ndimage.label(restricted_region_grids_with_defects)
    for i in range(1, nfeat_all + 1):
        defect_size = np.count_nonzero(label_all == i)
        if defect_size >= defect_size_cutoff:
            all_defect_sizes_cutoff.append(defect_size)

    total_grids = len(defect_data)

    defect_grids = np.sum(all_defect_sizes_cutoff)
    # print("sum of all defect:", defect_grids, "number of defects:", defect_grids_length)

    percent_defect_coverage = (defect_grids / total_grids) * 100

    # return avg defect size
    defect_avg_size = np.mean(all_defect_sizes_cutoff)

    # return total_num_defects_restricted_X, total_grids_restricted_X
    return defect_grids, total_grids, percent_defect_coverage, defect_avg_size


def surface_coverage_defects_cutoff(
    restriction, defect_size_cutoff, defect_cut_off, composition, configuration
):
    """
    Args:
        restriction: int of restricted area +/- from center of membrane that the defect coverage will be determined for.
        defect_size_cutoff: int the minimum area of a defect to be considered.
        defect_cut_off: int that defines the depth of a cutoff. This is a parameter when running the PackMem code.
        lipid_comp: str lipid composition of system.
        configuration: str of configuration of system.

    Returns:
        percent_coverage_defects: list containing avg defect coverage over frames, and std
        avg_defect_size_over_frames: list containing avg defect size overframes, and std
        n_frames: number of frames in the trajectory


    """

    # pulls the necessary data
    if configuration == "flat":
        time = "100ps-1800ext"
        analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}-equil-ext"
        analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
        folder_system = (
            analysis_defect_path
            / f"{composition}-{configuration}-production-stripped-ext1800-cutoff"
        )
    if configuration == "strain.2":
        time = "100ps-1800ext"
        analysis_path = f"/scratch/local/casakurai/asyn-phase-binding-data/analysis/defect-data-{time}-equil-ext"
        analysis_defect_path = Path(f"{analysis_path}/defect-cut-off-{defect_cut_off}A")
        folder_system = (
            analysis_defect_path
            / f"{composition}-{configuration}-production-stripped-centered-ext1800-cutoff"
        )

    defect_counts = np.load(folder_system / "output-defects-upper.npy")
    quad_centroids_np = np.load(folder_system / "output-centroids.npy")
    shape = np.load(folder_system / "output-shape.npy")
    X = np.load(folder_system / "output-Xshape.npy")
    Y = np.load(folder_system / "output-Yshape.npy")

    n_frames = len(quad_centroids_np)

    avg_defect_size_over_frames = {}

    total_defects_restricted_region_per_frame = np.empty(n_frames)
    total_grids_restricted_region_per_frame = np.empty(n_frames)
    percentage_defect_coverage_per_frame = np.empty(n_frames)
    avg_defect_size_per_frame = np.empty(n_frames)

    for frame in np.arange(0, n_frames):

        defect_counts_frame = defect_counts[frame]

        centroids_frame = quad_centroids_np[frame, :, 0:1]

        masked_all_defects_num = np.where(defect_counts_frame < 1, 1, 0)

        total_per_frame, total_grids, percent_defect_coverage, avg_defect_size = (
            total_num_defect_restricted_X_cutoff(
                masked_all_defects_num,
                shape,
                centroids_frame,
                restriction,
                defect_size_cutoff,
            )
        )

        total_defects_restricted_region_per_frame[frame] = total_per_frame

        total_grids_restricted_region_per_frame[frame] = total_grids

        percentage_defect_coverage_per_frame[frame] = percent_defect_coverage

        # handles frames with no defects
        if np.isnan(avg_defect_size):
            avg_defect_size_per_frame[frame] = 0
        else:
            avg_defect_size_per_frame[frame] = avg_defect_size

    average_defect_coverage = np.mean(percentage_defect_coverage_per_frame)

    # print(f"{composition} avg defect coverage {average_defect_coverage}")
    std_defect_coverage = np.std(percentage_defect_coverage_per_frame)

    average_defect_size = np.mean(avg_defect_size_per_frame)
    average_defect_size_std = np.std(avg_defect_size_per_frame)

    shape = percentage_defect_coverage_per_frame.shape
    num_frames = shape[0]

    percent_coverage_defects = [
        average_defect_coverage,
        std_defect_coverage,
        num_frames,
    ]
    avg_defect_size_over_frames = [
        average_defect_size,
        average_defect_size_std,
        num_frames,
    ]
    percent_coverage_defects_per_frame = percentage_defect_coverage_per_frame

    return (
        percent_coverage_defects,
        avg_defect_size_over_frames,
        percent_coverage_defects_per_frame,
        n_frames,
    )


def frame_extraction_avg_surface_defect_coverage(
    defect_restricted_avg_system, n_frames, percent_coverage_defects_per_frame
):
    """
    Compiles a list of frames that is +/- .01 away from the avg defect coverage over the entire membrane

    Args:
    defect_restricted_avg_system: list [avgerage defect size, std]
    n_frames: int number of frames in the trajectory
    percent_coverage_defects_per_frame: list of average defect coverage per frame in the trajectory

    Returns:
    max_frames_avg_defect_coverage_per_system (int): A single frame that represents the avg coverage of defects across the entire membrane for each lipid composition.

    """

    list_frames_close_avg_defect = []
    avg_defect = defect_restricted_avg_system[0]
    upper_bound = avg_defect + 0.01
    lower_bound = avg_defect - 0.01
    for frame in np.arange(0, n_frames):
        defect_coverage = percent_coverage_defects_per_frame[frame]
        if (
            lower_bound < defect_coverage < upper_bound
        ):  # storing frames close to the average defect coverage for entire membrane
            list_frames_close_avg_defect.append(frame)

    configurations_types_frames_to_avg_defect = list_frames_close_avg_defect

    max_frames_avg_defect_coverage_per_system  = max(configurations_types_frames_to_avg_defect)
    return max_frames_avg_defect_coverage_per_system
