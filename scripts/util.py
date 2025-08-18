
from pathlib import Path


#list
simulations = [
    1,
    2,
]

#dictionary
system_compositions = {
    1: {
        "DOPC": 100,
    },
    2: {
        "DPPC": 100,
    }
}


#Paths
base_path = Path("/scratch/casakurai/asyn-phase-binding-data/systems")

sim_path = Path("/scratch/casakurai/asyn-phase-binding-data/systems")
 
mdp_path = ("/scratch/casakurai/asyn-phase-binding-data/systems/mdps")

analysis_path = Path("/scratch/casakurai/asyn-phase-binding-data/analysis")

figures_path = Path("/scratch/casakurai/asyn-phase-binding-data/scripts/figures")

membrane_sel = "resname DPPC DOPC"
po4_sel = " name PO4 GL1 GL2 "
po4_only = "name PO4"
acyl_sel = "name C1A D2A C3A C4A C1B D2B C3B C4B"