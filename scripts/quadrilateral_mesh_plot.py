import pickle
import numpy as np
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter
from util import analysis_path

path = analysis_path / "curvature_selection" / "mesh"

# Load the data
with open(path / "membrane_curvature_smooth.pickle", "rb") as handle:
    mc = pickle.load(handle)

# Extract z surfaces
z_data_upper = mc.results["z_surface"]["upper"]
z_data_lower = mc.results["z_surface"]["lower"]

# Compute bin edges from x_range, y_range, and number of bins
x_edges = np.linspace(mc.x_range[0], mc.x_range[1], mc.n_x_bins + 1)
y_edges = np.linspace(mc.y_range[0], mc.y_range[1], mc.n_y_bins + 1)

# Compute bin centers
x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
X, Y = np.meshgrid(x_centers, y_centers, indexing="ij")

frame = 0  # choose frame
z_upper = z_data_upper[frame]
z_lower = z_data_lower[frame]

# # Optional: smooth the surfaces for better visualization
# z_upper_smooth = gaussian_filter(z_upper, sigma=1)
# z_lower_smooth = gaussian_filter(z_lower, sigma=1)

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
