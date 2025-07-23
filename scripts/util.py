
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

membrane_sel = "resname DPPC DOPC"
po4_sel = " name PO41 PO42 GLC "