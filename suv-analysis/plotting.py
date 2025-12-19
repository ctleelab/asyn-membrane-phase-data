import numpy as np
import pandas as pd
import ctleelab_plothelper.plothelpers as ph
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import matplotlib.cm as cm
from matplotlib.colors import TwoSlopeNorm

import pdb


df = pd.read_csv("curvatures.csv")


def plot_curve(df, i):
    curve = df.loc[df["Curve Name"] == f"CURVE {i}"]
    data = curve[
        [
            "X-Coordinate (um)",
            "Y-Coordinate (um)",
            "Point Curvature (um-1)",
            "Point Curvature Sign",
        ]
    ].to_numpy()

    # print(data)

    x = data[:, 0]
    y = data[:, 1]
    x -= np.mean(x)
    y -= np.mean(y)

    r = np.sqrt(x**2 + y**2)
    theta = np.mod(np.rad2deg(np.arctan2(y, x)), 360)
    curvature = data[:, 2] * data[:, 3]  # micron^-1
    # curvature_radius = np.divide(1, curvature) * 1000  # nm
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.true_divide(1, curvature)
        c[c == np.inf] = 999999
    curvature_radius = c * 1000  # nm

    cdf = pd.DataFrame(
        {"theta": theta, "curvature": curvature, "curvature_radius": curvature_radius}
    )
    cdf.sort_values("theta", inplace=True)

    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):

        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=1.5)

        ax.plot(cdf["theta"], cdf["curvature"])
        ax.set_xlabel("Theta ($\circ$)")
        ax.xaxis.set_major_locator(MultipleLocator(60))
        ax.set_ylabel("Curvature (μm$^{-1}$)")
        fig.savefig(f"curve_{i}.png")

        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=1.5)

        # ax.scatter(cdf["theta"], np.abs(cdf["curvature_radius"]), marker=".")
        ax.plot(cdf["theta"], np.abs(cdf["curvature_radius"]))
        ax.set_xlabel("Theta ($\circ$)")
        ax.xaxis.set_major_locator(MultipleLocator(60))
        ax.set_ylabel("Radius of Curvature (nm)")
        ax.set_ylim(0, 1000)

        fig.savefig(f"curve_{i}_radius.png")

        # Geometry plot: map curvature to colorbar on (x, y)
        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=1.5, subheight=1.5,   rmargin_scale = 2)

        norm = TwoSlopeNorm(
            vmin=np.nanmin(curvature), vcenter=0.0, vmax=np.nanmax(curvature)
        )

        sc = ax.scatter(x, y, c=curvature, cmap=cm.coolwarm, s=10, norm=norm)

        cbar = ph.add_fixed_colorbar(sc, ax=ax, aspect=20, pad=0.05)
        cbar.set_label("Curvature (μm$^{-1}$)")

        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")
        ax.set_aspect("equal", adjustable="datalim")
        fig.savefig(f"geom_{i}.png")


for i in range(1, 5):
    plot_curve(df, i)
