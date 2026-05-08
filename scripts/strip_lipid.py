import util as util
import MDAnalysis as mda 
import subprocess 
from util import analysis_path

strain = "flat"
size ="large"

systems = {
    "1":"DOPC",
    "2": "DPPC",
    "3": "DOPC-DOPS",
    "4": "DPPC-DPPS",
    "5": "DPPC-DPPA",
    "6": "DOPC-DOPA"
}


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

for system, lipid in systems.items(): 
    
    file_gro_path = analysis_path/strain/f"system{system}"/f"prod"
    input_file_xtc = file_gro_path /"production7.9.xtc"
    input_file_gro = file_gro_path /"production7.9.gro"
    output_file_name = file_gro_path / f"{lipid}-{strain}-production7.9-stripped"

    _strip_trajectory(input_file_gro, input_file_xtc)