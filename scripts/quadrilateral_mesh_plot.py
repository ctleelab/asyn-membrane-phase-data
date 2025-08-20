import pickle
import numpy as np
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter
from util import analysis_path

path = analysis_path / "curvature_selection" / "mesh"

# Load the data
with open(path / "membrane_curvature_smooth.pickle", "rb") as handle:
    mc = pickle.load(handle)

frame = 0  # choose frame
z_upper = mc["Z_upper"][frame]
z_lower = mc["Z_lower"][frame]
X = mc["X"]
Y = mc["Y"]
print(X.shape, Y.shape, z_upper.shape, z_lower.shape)
nx, ny = z_upper.shape
print("Number of bins in X:", nx)
print("Number of bins in Y:", ny)

# Create Plotly figure
fig = go.Figure()

#set min and max height (z-cord) between the leaflets and use that to set the color scale 
print("lowermin", z_lower.min())
print("uppermax", z_upper.max())

z_min = z_lower.min()
z_max = z_upper.max()


fig.add_trace(
    go.Surface(
        z=z_upper,
        x=X,
        y=Y,
        colorscale="Viridis",
        cmin = z_min,
        cmax = z_max,
        name="Upper leaflet",
        opacity=0.8,
    )
)

fig.add_trace(
    go.Surface(
        z=z_lower,
        x=X,
        y=Y,
        colorscale="Viridis",
        cmin = z_min,
        cmax = z_max,
        name="Lower leaflet",
        opacity=0.8,
    )
)

fig.update_layout(
    title=f"Membrane Leaflets Frame {frame}",
    scene=dict(
        xaxis_title="X",
        yaxis_title="Y",
        zaxis_title="Z height",
        aspectmode="data",  # keep axes aspect ratio consistent
    ),
)

# Save as standalone HTML
fig.write_html("membrane_leaflets_frame2501_nogaussian.html")
print("Interactive HTML saved as 'membrane_leaflets_frame2501_nogaussian.html'")
