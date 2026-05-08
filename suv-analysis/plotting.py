import numpy as np
import pandas as pd
import ctleelab_plothelper.plothelpers as ph
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import matplotlib.cm as cm
from matplotlib.colors import TwoSlopeNorm

import pdb
import csv
from scipy import stats


lipids = ["DOPC", "DPPC", "DOPC_adjusted", "DPPC_adjusted", "DPPC_DPPS", "DOPC_DOPS"]


def polygon_area(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate the signed area of a polygon given its vertices using the Shoelace formula.
    """
    return 0.5 * (np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def polygon_perimeter(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate the perimeter of a closed polygon given its vertices.
    """
    dx = x - np.roll(x, 1)
    dy = y - np.roll(y, 1)
    return np.sum(np.sqrt(dx**2 + dy**2))


def plot_curve(df, i, lipid, format="png"):
    curve = df.loc[df["Curve Name"] == f"CURVE {i}"]

    if len(curve) == 0:
        return

    data = curve[
        [
            "X-Coordinate (um)",
            "Y-Coordinate (um)",
            "Point Curvature (um-1)",
            "Point Curvature Sign",
            "Curve Length (um)",
        ]
    ].to_numpy()

    x = data[:, 0] * 1000
    y = data[:, 1] * 1000
    x -= np.mean(x)
    y -= np.mean(y)

    area = polygon_area(x, y)
    print(f"{lipid} curve {i} area: {area} nm^2")

    ideal_radius = np.sqrt(np.abs(area) / np.pi) / 1000  # micron
    ideal_curvature = 1 / ideal_radius
    print(f"Ideal curvature for {lipid} curve {i}: {ideal_curvature} micron^-1")

    # perimeter = polygon_perimeter(x, y) # length from polygon
    perimeter = data[0, 4] * 1000  # length from bezier curve nm
    # print(f"{lipid} curve {i} perimeter: {perimeter} nm, {data[0,4] * 1000} nm")

    r = np.sqrt(x**2 + y**2)
    theta = np.mod(np.rad2deg(np.arctan2(y, x)), 360)
    curvature = data[:, 2] * data[:, 3]  # micron^-1
    if area < 0:
        curvature = -curvature
        area = -area
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
        ax.set_xlabel(r"Theta ($\circ$)")
        ax.xaxis.set_major_locator(MultipleLocator(60))
        ax.set_ylabel("Curvature (μm$^{-1}$)")
        ax.set_title(f"{lipid} curve{i}")
        fig.savefig(f"{lipid}_{i}.{format}")
        plt.close(fig)

        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=1.5)

        # ax.scatter(cdf["theta"], np.abs(cdf["curvature_radius"]), marker=".")
        ax.plot(cdf["theta"], np.abs(cdf["curvature_radius"]))
        ax.set_xlabel(r"Theta ($\circ$)")
        ax.xaxis.set_major_locator(MultipleLocator(60))
        ax.set_ylabel("Radius of Curvature (nm)")
        ax.set_ylim(0, 1000)
        ax.set_title(f"{lipid} curve{i}")

        fig.savefig(f"{lipid}_{i}_radius.{format}")
        plt.close(fig)

        # Geometry plot: map curvature to colorbar on (x, y)
        fig, ax = ph.fixed_size_subplots(
            1, 1, subwidth=1.5, subheight=1.5, rmargin_scale=2
        )

        print(
            lipid,
            "curve",
            i,
            np.nanmax(curvature),
            np.nanmin(curvature),
            np.mean(curvature),
            np.std(curvature),
        )

        norm = TwoSlopeNorm(vmin=-300, vcenter=0.0, vmax=300)

        sc = ax.scatter(x, y, c=curvature, cmap=cm.PRGn, s=4, norm=norm)

        cbar = ph.add_fixed_colorbar(sc, ax=ax, aspect=20, pad=0.05, extend="both")
        cbar.set_label("Curvature (μm$^{-1}$)")
        cbar.locator = MultipleLocator(100)
        cbar.update_ticks()

        ax.set_xlabel("X (nm)")
        ax.set_ylabel("Y (nm)")

        ax.set_xlim(-35, 35)
        ax.set_ylim(-35, 35)

        ax.set_title(f"{lipid} curve{i}")
        ax.set_aspect("equal", adjustable="datalim")
        fig.savefig(f"{lipid}_{i}_geom.{format}")
        plt.close(fig)

        delta_ideal_curvature = curvature - ideal_curvature

        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=1.5)
        ax.plot(cdf["theta"], delta_ideal_curvature)
        ax.set_xlabel(r"Theta ($\circ$)")
        ax.xaxis.set_major_locator(MultipleLocator(60))
        ax.set_ylabel("Curvature Deviation (μm$^{-1}$)")
        ax.set_title(f"{lipid} curve{i} curvature deviation")
        fig.savefig(f"{lipid}_{i}_delta_ideal.{format}")
        plt.close(fig)

        return 4 * np.pi * area / (perimeter**2)

        delta_ideal_curvature = curvature - ideal_curvature

for lipid in lipids:
    df = pd.read_csv(f"traces_{lipid}.csv")
    for i in range(1, 6):
        polsby_popper_score = plot_curve(df, i, lipid)
        print(polsby_popper_score)
        print()
