import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ctleelab_plothelper.plothelpers as ph
from pathlib import Path



# load csv
time = "100ps-1800ext"
data_path = Path(f"/scratch/local/casakurai/asyn-phase-binding-data/Figures/defect-data-{time}-equil-ext")
restriction = 4
circle_fit_path = data_path/f"circle-fit"
curvature_data = circle_fit_path/f"all_systems_curvature_ext1800_{restriction}nmregion.csv"
df = pd.read_csv(curvature_data)

systems = df["system"]
x = np.arange(len(systems))

print(df)

############################
# Radius plot
############################

fig, ax = ph.fixed_size_subplots(1,1, subwidth=3, subheight=1.25)

ax.bar(
    x,
    df["avg_radii(nm)"],
    yerr=df["std_radii"],
    capsize=5
)

ax.set_xticks(x)
ax.set_xticklabels(systems, rotation=45)
ax.set_ylabel("Radius (nm)")
ax.set_title("Average Buckle Radius")


fig.savefig(f"{circle_fit_path}/radius_barplot_{restriction}nmregion.pdf",bbox_inches="tight")

############################
# Curvature plot
############################

fig, ax = ph.fixed_size_subplots(1,1, subwidth=3, subheight=1.25)

ax.bar(
    x,
    df["avg_curvature(um^-1)"],
    yerr=df["std_curvature"],
    capsize=5, 
    color = "#0000ff"
)

ax.set_xticks(x)
ax.set_xticklabels(systems, rotation=45)
ax.set_ylabel("κ (µm$^{-1}$)")
ax.set_title("Average Buckle Curvature Entire membrane")

fig.tight_layout()
fig.savefig(f"{circle_fit_path}/curvature_barplot{restriction}nmregion.pdf",bbox_inches="tight")
