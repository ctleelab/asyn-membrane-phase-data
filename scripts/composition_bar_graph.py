import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
import matplotlib.colors as mc
import colorsys

# Your compositions dictionary
compositions = { 
    1: {"DOPC":100, "DPPC":0, "DOPS":0, "DPPS":0, "DOPA":0, "DPPA":0},
    2: {"DOPC":0, "DPPC":100, "DOPS":0, "DPPS":0, "DOPA":0, "DPPA":0},
    3: {"DOPC":75, "DPPC":0, "DOPS":25, "DPPS":0, "DOPA":0, "DPPA":0},
    4: {"DOPC":0, "DPPC":75, "DOPS":0, "DPPS":25, "DOPA":0, "DPPA":0},
    5: {"DOPC":0, "DPPC":75, "DOPS":0, "DPPS":0, "DOPA":0, "DPPA":25},
    6: {"DOPC":75, "DPPC":0, "DOPS":0, "DPPS":0, "DOPA":25, "DPPA":0}
}

# Lipid order
lipid_names = ["DOPC", "DPPC", "DOPS", "DPPS", "DOPA", "DPPA"]

lipid_hatches = {
    "DOPC": None,
    "DPPC": None,   
    "DOPS": "////", 
    "DPPS": "////",
    "DOPA": "....",
    "DPPA": "....",
}

lipid_colors = {
    "DOPC": "#aaa9a9ff",   
    "DPPC": "#d2826cff",   
    "DOPS": "#aaa9a9ff",  
    "DPPS": "#d2826cff",   
    "DOPA": "#aaa9a9ff",   
    "DPPA":"#d2826cff",   
}

# Color helper
def lighten_color(color, amount=0.7):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

# Output folder
curr_fig_path = Path("Figures")
curr_fig_path.mkdir(exist_ok=True)

# Color palette
pal = sns.color_palette("colorblind", len(lipid_names))

# Plot styles
plot_styles = [
    ("ctleelab_plothelper.base", ""),
    ("ctleelab_plothelper.dark", "_dark"),
]

for style, style_ext in plot_styles:

    with plt.style.context(style):

        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=3)

        # Plot stacked bars
        for sim, composition in compositions.items():

            bottom = 0

            for i, lipid in enumerate(lipid_names):

                value = composition[lipid]

                if value == 0:
                    continue

                ax.bar(
                    sim,
                    value,
                    bottom=bottom,
                    color=lighten_color(lipid_colors[lipid], 0.8),
                    edgecolor="black",
                    linewidth=1,
                    width=0.8,
                    hatch=lipid_hatches[lipid],   # ← ADD THIS LINE
                    label=lipid if sim == 1 else ""
                )
                bottom += value

        # Formatting
        ax.set_xlabel("System")
        ax.set_ylabel("Composition (%)")

        ax.set_ylim(0, 100)
        ax.set_xticks(list(compositions.keys()))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, 1.2),
            ncols=len(lipid_names),
            frameon=False
        )

        fig.tight_layout()

        # Save
        fig.savefig(curr_fig_path / f"lipid_compositions{style_ext}.pdf")

        plt.show()
        plt.close(fig)
