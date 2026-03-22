# script not set up to only do one dimension type
# import libraries
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import MDAnalysis as mda

# turns the edr file into a numpy array


systems = ["5", "6"]
target_strain = 0.2
# compare strain using Lo from flat membranes vs Lo of buckled
system_strain_comparison = {}
compression_data_analysis = Path(
    "/scratch/local/casakurai/asyn-phase-binding-data/compression-data/analysis"
)


def compression_plot_set_Lo(system):
    # first grab the length of Lo (which is the average Lx for flat membrane in NPT)
    flat_large_path = Path(
        f"/scratch/local/casakurai/asyn-phase-binding-data/systems/system{system}-8x8x25/equil-large/"
    )
    flat_edr_file = flat_large_path / "equilibration6.7.edr"
    aux_flat = mda.auxiliary.EDR.EDRReader(flat_edr_file)

    terms_flat = aux_flat.get_data(["Box-X", "Box-Z", "Time"])

    flat_box_x_avg = np.mean(terms_flat["Box-X"])

    # then grab the starting Lo of x-dim of the compression simulation
    compression_large_path = Path(
        f"/scratch/local/casakurai/asyn-phase-binding-data/compression-data/system{system}"
    )
    compression_edr_file = compression_large_path / "large-compression.edr"
    aux_compression = mda.auxiliary.EDR.EDRReader(compression_edr_file)
    terms_compression = aux_compression.get_data(["Box-X", "Box-Z", "Time"])

    # determining the minimum Lo for compression simulation
    time0_idx = np.argmin(np.abs(terms_compression["Time"] - 0))
    print(time0_idx)
    compression_initial_dim_x = terms_compression["Box-X"][time0_idx]

    # from buckled
    compression_strain = (
        compression_initial_dim_x - terms_compression["Box-X"]
    ) / compression_initial_dim_x
    Time_1e6 = terms_compression["Time"] / 1000000

    # combine into dataframe
    df = pd.DataFrame(
        {
            "Time_1e6": Time_1e6,
            "Time": terms_compression["Time"],
            "Box-X": terms_compression["Box-X"],
            "Compression strain": compression_strain,
        }
    )

    # find the index where the compression strain is closest to the target strain
    # this might not be entirely accurate
    strain_idx = np.argmin(np.abs(df["Compression strain"] - target_strain))

    time_at_target_strain = df["Time"].iloc[strain_idx]
    compression_at_target_strain = df["Compression strain"].iloc[strain_idx]
    Lt_at_target_strain = df["Box-X"].iloc[strain_idx]

    # new strain based on Lo of flat NPT
    compression_strain_adj = (
        flat_box_x_avg - terms_compression["Box-X"]
    ) / flat_box_x_avg
    df_new = pd.DataFrame(
        {
            "Time_1e6": Time_1e6,
            "Time": terms_compression["Time"],
            "Box-X": terms_compression["Box-X"],
            "Compression strain": compression_strain,
            "adjusted compression strain": compression_strain_adj,
        }
    )
    adjusted_target_strain = df_new["adjusted compression strain"].iloc[strain_idx]

    # correct target strain
    idx_at_targert_strain_new = np.argmin(
        np.abs(df_new["adjusted compression strain"] - target_strain)
    )
    value_idx_at_target_strain_new = df_new["adjusted compression strain"].iloc[
        idx_at_targert_strain_new
    ]
    time_at_target_strain_new = df_new["Time"].iloc[idx_at_targert_strain_new]

    print(f"system{system}, Lo based on flat NVT is {adjusted_target_strain}")
    differenece_between_strains = abs(target_strain - adjusted_target_strain)
    print(
        f"The adjusted strain based on Lo of avg Lx of flat NVT would be {adjusted_target_strain}"
    )

    system_strain_comparison[system] = [
        target_strain,
        adjusted_target_strain,
        differenece_between_strains,
        compression_at_target_strain,
    ]

    return df_new, value_idx_at_target_strain_new, time_at_target_strain_new


# plot the differences between the strains
def bar_plot_strains(data):
    systems = [
        "DOPC",
        "DPPC",
        "DOPC-DOPS",
        "DPPC-DPPS",
        "DPPC-DPPA",
        "DOPC-DOPA",
    ]
    differences = []
    new_strains = []
    actual_strains = []
    for key, value in data.items():
        difference_between_strains = value[2]
        new_strain = value[1]
        actual_strain = value[3]
        print(actual_strain)
        differences.append(difference_between_strains)
        new_strains.append(new_strain)
        actual_strains.append(actual_strain)

    plt.bar(systems, differences)

    plt.xlabel("System")
    plt.ylabel("Difference Between Strains")
    plt.title("Strain Difference by System")

    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(compression_data_analysis / "strain_comparison.png")
    plt.close()

    plt.bar(systems, new_strains)
    plt.xlabel("System")
    plt.ylabel("Correct strains")
    plt.title("Strains based on Lo of flat NPT avg Lx")
    plt.axhline(y=0.2, color="red", linestyle="--")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(compression_data_analysis / "adjusted_strains.png")
    plt.close()

    plt.bar(systems, actual_strains)
    plt.xlabel("System")
    plt.ylabel("actual strains")
    plt.title("Actual strains for those closest to strain .2")
    plt.axhline(y=0.2, color="red", linestyle="--")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(compression_data_analysis / "actual_strain_at_strain2.png")
    plt.close()


# bar_plot_strains(system_strain_comparison)


def compression_plot_strain(
    compression_data, system, value_idx_at_targert_strain_new, time_at_target_strain_new
):
    # combined plot time vs compression
    plt.plot(compression_data["Time_1e6"], compression_data["Compression strain"])
    plt.xlabel("time (ps)")
    plt.ylabel("compression strain")
    plt.title(
        f"system{system} target-strain{target_strain} real-strain {value_idx_at_targert_strain_new:.4g} at {time_at_target_strain_new}ps "
    )

    # Set major ticks every 0.1 on y-axis
    plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))

    # Draw horizontal red line at compression = 0.1
    plt.axhline(y=0.1, color="red", linestyle="--", linewidth=1)
    plt.axhline(y=0.2, color="orange", linestyle="--", linewidth=1)
    plt.axhline(y=0.6, color="green", linestyle="--", linewidth=1)

    plt.tight_layout()
    plt.savefig(
        compression_data_analysis
        / f"compression_{target_strain}strain_system{system}_{time_at_target_strain_new}ps.png"
    )
    plt.close()


for system in systems:
    compression_data, value_idx_at_target_strain_new, time_at_target_strain_new = (
        compression_plot_set_Lo(system)
    )
    compression_plot_strain(
        compression_data,
        system,
        value_idx_at_target_strain_new,
        time_at_target_strain_new,
    )
