
from pathlib import Path


#list
simulations = [
    1,
    2,
    3,
    4,
    5,
    6
]

#dictionary
system_compositions = {
    1: {
        "DOPC": 100,
    },
    2: {
        "DPPC": 100,
    },
    3: {
        "DOPC": 75,
        "DOPS":25
    },
    4: {
        "DPPC":75,
        "DPPS": 25
    },
    5: {
        "DPPC": 75,
        "DPPA": 25
    },
    6: {
        "DOPC": 75,
        "DOPA": 25
    }
}


#Paths
base_path = Path("/scratch/casakurai/asyn-phase-binding-data/systems")

sim_path = Path("/scratch/casakurai/asyn-phase-binding-data/systems")
 
mdp_path = ("/scratch/casakurai/asyn-phase-binding-data/systems/mdps")

analysis_path = Path("/scratch/casakurai/asyn-phase-binding-data/simulations")

figures_path = Path("/scratch/casakurai/asyn-phase-binding-data/scripts/figures")

membrane_sel = "resname DOPC DPPC DOPS DPPS DPPA DOPA"
po4_sel = " name PO4 GL1 GL2 "
po4_only = "name PO4"
acyl_sel = "name C1A D2A C3A C4A C1B D2B C3B C4B"