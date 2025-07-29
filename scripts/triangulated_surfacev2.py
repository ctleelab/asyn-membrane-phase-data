#cluster lipids by leaflets (Leaflet finder)
#create a function that represents the cloud of point(PO4 beads)
#then create vectors for each of the cubic/linear functions
#then take the Laplacian of the vector field which will represent the surface 

import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder

import pathlib as Path 
import util
from util import analysis_path 
from util import figures_path
import numpy as np 
import matplotlib.pyplot as plt

#paths that are unchanged
curvature_selection_path = analysis_path/f"curvature_selection"
figures_triangualted_surface = figures_path/"triangulated_surface"
figures_triangualted_surface.mkdir(parents=True, exist_ok=True)


lipids = ['DOPC']



def leaflet_def(lipid: str) :

    u = mda.Universe(f"{system_path}/large-system.top", f"{equil_path}/equilibration6.6.gro", topology_format = "ITP")

    #defining the upper and lower leaflet with PO4 beads 
    L = mda.analysis.leaflet.LeafletFinder(u,'name PO4')
    leaflet0 = L.groups(0)
    leaflet1 = L.groups(1)

    leaflets = [leaflet0, leaflet1]
    return(leaflets)

#get the normal vector for the PO4 bead that aligns with the lipid tail (how to define the lipid tail?)
#there's two acyl tails..
#need to complete for each residue, and complete for each leaflet 

def COM_def(residue:str):
    acyl_chain = ['C1A','D2A', 'C3A' ,'C4A' ,'C1B' ,'D2B' ,'C3B' ,'C4B']
    COM_atoms_df = []
    total_mass = 0
    for atom in residue.atoms:
        #center of mass of each atom in acyl chain
        if atom.name in acyl_chain:
            COM_atom = atom.position * atom.mass
            COM_atoms_df.append(COM_atom)
    total_mass = sum([atom.mass for atom in residue.atoms])
    COM_residue = sum(COM_atoms_df)/total_mass

    return COM_residue



#distinguishes the two leaflets
for lipid in lipids:
    system_path = curvature_selection_path/f"system1-8x8x25-{lipid}"
    equil_path = system_path/f"equil"
    leaflet_def(lipid)


#creating a normal vector between PO4 and COM of acyl chain 
normal_vectors = []
COM_residues_df = []
PO4_position = []
leaflets = leaflet_def(lipids)
for leaflet in leaflets:
    for residue in leaflet.residues:
        PO4_atom = residue.atoms.select_atoms(util.po4_only)
        if len(PO4_atom) == 1:
            PO4_pos = PO4_atom[0].position
            PO4_position.append(PO4_pos)
        else: 
            continue 
        if not PO4_atom:
            continue

        COM_tail = COM_def(residue)

        vec = COM_tail - PO4_pos
        #computes the magnitude of the vector
        norm = np.linalg.norm(vec)
        #divide by the magnitude so it's a unit vector of 1
        vec /= norm
        vec = np.array(vec)
        normal_vectors.append(vec)

#representing vectors on a graph
PO4_position = np.array(PO4_position)  # shape (N,3)
vectors = np.array(normal_vectors)        # shape (N,3), normalized


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.quiver(
    PO4_position[:,0], PO4_position[:,1], PO4_position[:,2],
    vectors[:,0], vectors[:,1], vectors[:,2],
    length=3.0, normalize=True, color='r'
)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title("Vectors Originating at PO4 Positions")
plt.savefig(figures_triangualted_surface/f"DOPC-equilibration6.6.png")
plt.show()





    

