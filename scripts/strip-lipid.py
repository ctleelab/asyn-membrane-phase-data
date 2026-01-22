import util as util
import MDAnalysis as mda 
import subprocess 
from util import analysis_path

system = 2
shape = "buckled"
lipid = "DPPC"
strain = "0.15strain"
file_gro_path = analysis_path/"curvature_selection"/"NVT"/shape/f"system{system}-8x8x25-{lipid}-{strain}-NVT"/"prod"
input_file_xtc = file_gro_path /"production7.7.xtc"
input_file_gro = file_gro_path /"production7.7.gro"
output_file_name = file_gro_path / "production7.7_stripped"

def _strip_trajectory(gro,xtc) -> None:

    # Input .gro (already provided)

    # Load your system
    u = mda.Universe(gro,xtc)

    # Select only lipid residues (replace with your lipid names)
    lipids = u.select_atoms("resname DPPC DOPC DOPS DPPS DOPA DPPA")  

    # Write to a new .gro file
    lipids.write(f"{output_file_name}.gro")

    # Write a stripped trajectory (coordinates for only lipids)
    with mda.Writer(f"{output_file_name}.xtc", lipids.n_atoms) as W:
        for ts in u.trajectory:
            W.write(lipids)


_strip_trajectory(input_file_gro, input_file_xtc)