import MDAnalysis as mda
from pathlib import Path
import numpy as np
from util import analysis_path

shape = "buckled"
strain = "0.2strain"


#inputs
systems = { "1": "DOPC",
           "2": "DPPC",
    "3": "DOPC-DOPS",
    "4": "DPPC-DPPS",
    "5": "DPPC-DPPA",
    "6": "DOPC-DOPA"
}


def centered_membrane(gro,xtc):
    u = mda.Universe(gro, xtc)
    #select lipids 
    lipids = u.select_atoms("resname DOPC DPPC DOPS DPPS DOPA DPPA")
    all_atoms = u.select_atoms("all")


    # --- First frame for GRO ---
    ts = u.trajectory[0]  # first frame

    # Center first frame based on peak lipids
    z_positions = lipids.positions[:, 2]
    x_positions = lipids.positions[:, 0]
    peak_z_value = np.max(z_positions) #finds the max z position out of all lipids
    peak_z_value_idx = np.argmax(z_positions)
    peak_lipid_x = lipids.positions[peak_z_value_idx, 0]
    
    z_threshold = peak_z_value - 4 #4nm below the peak, want it to be ~thickness of the membrane
    x_threshold_min = peak_lipid_x - 1
    x_threshold_max = peak_lipid_x + 1
    
    #create two boolean mask for lipids that are within the z threshold and for the x threshold
    #choose lipids that are within 10A of the peak z position
    peak_lipids_mask_z = z_positions >= z_threshold 
    peak_lipid_mask_x =  (x_positions >= x_threshold_min) & (x_positions <= x_threshold_max)

    #combined mask 
    total_peak_lipids_mask = peak_lipids_mask_z & peak_lipid_mask_x

    #define the lipids that are within the specified peak regions
    peak_lipids = lipids[total_peak_lipids_mask]

    #calculate the center of mass of the peak lipids
    com_peak = peak_lipids.center_of_mass()

    #shift the box center based on the com of the peak lipids
    box_center = np.array([ts.dimensions[0]/2, ts.dimensions[1]/2, ts.dimensions[2]/2]) #center the com at (0,0,0)
    shift = box_center - com_peak
    all_atoms.positions = all_atoms.positions.copy() + shift

    #PBC
    all_atoms.wrap(compound="atoms")


    # --- Save first frame as GRO ---
    all_atoms.write(f"{output_file_name}.gro")


    with mda.Writer(f"{output_file_name}.xtc", all_atoms.n_atoms) as W:
        for ts in u.trajectory:
            z_positions = lipids.positions[:,2] #pulls out the z positions of all lipids
            x_positions = lipids.positions[:,0] #pulls out the x positions of all lipids
            peak_z_value = np.max(z_positions) #finds the max z position out of all lipids

            #finds the indx of the the lipid wiht the hightest z value
            peak_z_value_idx = np.argmax(z_positions)

            #pulls out the (x,y) of the max lipid
            peak_lipid_x = lipids.positions[peak_z_value_idx, 0]

            #select lipids within a cuttoff of the highest lipid's z-coord
            z_threshold = peak_z_value - 4 #4nm below the peak, want it to be ~thickness of the membrane

            #select lipids within a cutoff of the highest lipid's x-coord to avoid two peaks being selected
            # +/- 1nm in the x direction
            x_threshold_min = peak_lipid_x - 1
            x_threshold_max = peak_lipid_x + 1

            #create two boolean mask for lipids that are within the z threshold and for the x threshold
            #choose lipids that are within 10A of the peak z position
            peak_lipids_mask_z = z_positions >= z_threshold 
            peak_lipid_mask_x =  (x_positions >= x_threshold_min) & (x_positions <= x_threshold_max)

            #combined mask 
            total_peak_lipids_mask = peak_lipids_mask_z & peak_lipid_mask_x

            #define the lipids that are within the specified peak regions
            peak_lipids = lipids[total_peak_lipids_mask]

            #calculate the center of mass of the peak lipids
            com_peak = peak_lipids.center_of_mass()

            #shift the box center based on the com of the peak lipids
            box_center = np.array([ts.dimensions[0]/2, ts.dimensions[1]/2, ts.dimensions[2]/2]) #center the com at (0,0,0)
            shift = box_center - com_peak
            all_atoms.positions = all_atoms.positions.copy() + shift

            #PBC
            all_atoms.wrap(compound="atoms")

            #write out the frame
            W.write(all_atoms)


for system, lipid in systems.items(): 
    file_gro_path = analysis_path/"curvature_selection"/"NVT"/shape/f"system{system}-8x8x25-{lipid}-{strain}-NVT"/"prod"
    input_file_xtc = file_gro_path /"production7.8_stripped.xtc"
    input_file_gro = file_gro_path /"production7.8_stripped.gro"
    output_file_name = file_gro_path / "production7.8_stripped_centered"
    centered_membrane(input_file_gro, input_file_xtc)
