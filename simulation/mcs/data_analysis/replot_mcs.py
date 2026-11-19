"""Replot saved publication data: python replot_mcs.py plot_data/h2o/file.npz.
Place beside the updated mcs_helper.py. No ROOT input or model calculation needed.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

STYLE = {"font.size": 16, "axes.labelsize": 11, "legend.fontsize": 8,
         "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
         "lines.linewidth": 1.4, "pdf.fonttype": 42, "ps.fonttype": 42}


def plot_excess(data):
    fig, (ax, res) = plt.subplots(
        2, 1, sharex=True, figsize=(15 / 2.54, 14 / 2.54),
        gridspec_kw={"height_ratios": [3, 1]}, constrained_layout=True,
    )
    for key in data:
        if not key.endswith("_x_mm"):
            continue
        prefix = key[:-4]
        d = {k[len(prefix):]: v for k, v in data.items() if k.startswith(prefix)}
        single = prefix.startswith("single_")
        color, line, sampled = ("#0072B2", "-", "--") if single else ("#D55E00", "-.", ":")
        name = "Single Gaussian" if single else "Gaussian mixture"
        ax.plot(d["x_mm"], d["density"], color=color, linestyle=line, label=str(d["label"]))
        if "sample_density" in d:
            ax.stairs(d["sample_density"], d["sample_edges_mm"], baseline=None,
                      color=color, linestyle=sampled, linewidth=1.4, label=f"{name}, sampled")
        if "edges_mm" not in d:
            continue
        edges = d["edges_mm"]
        centers = (edges[:-1] + edges[1:]) / 2
        if bool(d["show_geant4"]):
            ax.errorbar(centers, d["g4_density"], yerr=d["g4_error"], fmt="o",
                        markersize=3, markerfacecolor="none", color="black",
                        elinewidth=0.8, capsize=0, label="Geant4")
            res.fill_between(centers, -d["g4_error"], d["g4_error"],
                             color="0.85", linewidth=0, zorder=0,
                             label=r"Geant4 stat. uncertainty ($\pm1\sigma$)")
        if "residual" in d:
            res.plot(centers, d["residual"], color=color, linestyle=line, linewidth=1.4)
        if "sample_residual" in d:
            res.stairs(d["sample_residual"], edges, baseline=None,
                       color=color, linestyle=sampled, linewidth=1.4)
    ax.set_ylabel(r"Probability density / mm$^{-1}$")
    res.set(xlabel=r"Excess path length $\Delta r$ / mm",
            ylabel="Residual\n" + r"/ mm$^{-1}$")
    ax.set_axisbelow(True)
    for panel in (ax, res):
        panel.grid(True, which="major", color="0.9", linestyle="-", linewidth=0.7)
    res.axhline(0, color="0.4", linewidth=0.8, zorder=0)
    ax.legend()
    if res.get_legend_handles_labels()[0]:
        res.legend(frameon=False, fontsize=8, loc="best")
    ax.set(xlim=data["xlim"], ylim=data["ylim"])
    res.set_ylim(data["res_ylim"])
    return fig


def replot(filename):
    with np.load(filename, allow_pickle=False) as archive:
        data = dict(archive)
    kind = str(data.pop("kind"))
    with plt.rc_context(STYLE):
        if kind == "excess":
            return {"excess": plot_excess(data)}
        if kind != "reach":
            raise ValueError(f"Unknown plot kind: {kind}")
        from mcs_helper import plot_reach_and_stopping
        data.pop("E0", None)
        reach, stop = plot_reach_and_stopping(**data)
        return {"reach": reach, "stopping": stop}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path)
    parser.add_argument("--save", action="store_true", help="Also save PDFs beside the NPZ")
    args = parser.parse_args()
    for name, fig in replot(args.file).items():
        if args.save:
            fig.savefig(args.file.with_name(f"{args.file.stem}_{name}.pdf"))
    plt.show()
