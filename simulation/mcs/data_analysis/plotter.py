"""Plot a selected layer's signed cumulative projected angle (input: radians).

Standalone:
    python plot_scattering_core.py h2oproj.root --layer 100

Inside the existing per-layer loop, AFTER CumSigmaCore[i] is assigned:
    if i == PLOT_LAYER_INDEX:
        plot_scattering_core(CumScatteringAngle[rows], i, z[i],
                             sigma_rad=CumSigmaCore[i])

The Gaussian is normalized to unit area and centred at zero, as in the
analytical model. No Gaussian fit, peak rescaling, or tail rejection is used.
"""
from pathlib import Path
from math import erfc, sqrt
import argparse
import numpy as np
import matplotlib.pyplot as plt


def plot_scattering_core(angles_rad, layer_index, depth_cm, sigma_rad=None,
                         output_dir="scattering_core", show=True):
    """Return (figure, axes); save SVG, PDF and 220-dpi PNG.

    sigma_rad=None reproduces the script's percentile width for >=100 hits.
    With fewer hits, pass the width actually used by the model explicitly:
    its previous-layer fallback cannot be inferred from this layer alone.
    Nonfinite entries are excluded and their count is reported.
    """
    raw = np.asarray(angles_rad, dtype=float)
    if raw.ndim != 1:
        raise ValueError("Expected a 1D array of signed projected angles in radians.")
    finite = raw[np.isfinite(raw)]
    dropped = raw.size - finite.size
    if finite.size == 0:
        raise ValueError("This layer has no finite angles.")
    if sigma_rad is None:
        if finite.size < 100:
            raise ValueError("Fewer than 100 finite hits: pass the model's sigma_rad explicitly.")
        low, high = np.percentile(finite, [15.8655, 84.1345])
        sigma_rad = 0.5 * (high - low)
    if not np.isfinite(sigma_rad) or sigma_rad <= 0:
        raise ValueError("The model core width must be finite and positive.")

    angles = finite * 1000.0  # radians -> mrad
    sigma = float(sigma_rad) * 1000.0
    n = angles.size
    tail_fraction = np.mean(np.abs(angles) > 3 * sigma)
    gaussian_tail = erfc(3 / sqrt(2))
    sample_std = np.std(angles)
    full_limit = max(8 * sigma, float(np.max(np.abs(angles))) * 1.015)
    blue, orange, ink = "#21618C", "#D55E00", "#243746"

    style = {"font.family": "DejaVu Sans", "font.size": 11,
             "axes.labelsize": 12, "axes.titlesize": 13,
             "axes.spines.top": False, "axes.spines.right": False,
             "axes.edgecolor": "#B2BEC6", "axes.labelcolor": ink,
             "text.color": ink, "xtick.color": ink, "ytick.color": ink,
             "svg.fonttype": "none", "pdf.fonttype": 42}
    with plt.rc_context(style):
        fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.7))
        fig.subplots_adjust(left=.075, right=.975, bottom=.27, top=.79, wspace=.27)
        fig.suptitle("Cumulative scattering: Gaussian core and tails",
                     x=.075, y=.965, ha="left", fontsize=19, fontweight="bold")
        fig.text(.075, .895,
                 f"Layer {layer_index}  |  Depth {depth_cm:.4f} cm  |  "
                 f"{n:,} finite hits  |  " + r"$\sigma_{\rm core}$" + f" = {sigma:.3g} mrad",
                 fontsize=11.5)

        for ax, limit, logarithmic in zip(axes, [4*sigma, full_limit], [False, True]):
            # Normalize by ALL finite hits, including those outside the core view.
            bins = int(np.clip(np.ceil(2*limit/(sigma/8)), 64, 2400))
            edges = np.linspace(-limit, limit, bins + 1)
            counts, _ = np.histogram(angles, bins=edges)
            density = counts / (n * np.diff(edges))
            values = np.where(counts > 0, density, np.nan) if logarithmic else density
            ax.stairs(values, edges, color=blue, linewidth=1.5,
                      label="Geant4 distribution", zorder=3)
            if not logarithmic:
                ax.stairs(density, edges, fill=True, color=blue, alpha=.10)
            x = np.linspace(-limit, limit, 4001)
            gaussian = np.exp(-.5*(x/sigma)**2) / (sqrt(2*np.pi)*sigma)
            ax.plot(x, gaussian, color=orange, linewidth=2.3,
                    label="Gaussian core model", zorder=4)
            ax.axvspan(-sigma, sigma, color=orange, alpha=.075, zorder=0)
            for threshold in [-3*sigma, 3*sigma]:
                ax.axvline(threshold, color=ink, linestyle=(0, (3, 4)),
                           linewidth=.9, alpha=.55)
            ax.set_xlim(-limit, limit)
            ax.set_xlabel(r"Cumulative projected angle $\theta_x$ / mrad")
            ax.set_ylabel(r"Probability density / mrad$^{-1}$")
            ax.grid(axis="y", alpha=.16)
            ax.set_axisbelow(True)
            if logarithmic:
                positive = density[density > 0]
                ax.set_yscale("log")
                ax.set_ylim(positive.min()*.35,
                            max(positive.max(), 1/(sqrt(2*np.pi)*sigma))*2)
                ax.set_title("Tails · logarithmic scale", loc="left", fontweight="bold", pad=12)
            else:
                ax.set_ylim(bottom=0)
                ax.set_title("Core · linear scale", loc="left", fontweight="bold", pad=12)
                ax.legend(frameon=False, fontsize=10, loc="upper left")

        fig.text(.075, .145,
                 r"Beyond $\pm3\sigma_{\rm core}$:  "
                 f"Geant4 {100*tail_fraction:.2f}%   |   Gaussian {100*gaussian_tail:.2f}%"
                 f"     ·     Full-sample std / core width = {sample_std/sigma:.2f}",
                 fontsize=11, fontweight="bold")
        fig.text(.075, .055,
                 "Gaussian centred at zero; width from the model's percentile-core prescription.\n"
                 "Shading: ±1 core width. Dashed lines: ±3 core widths (diagnostic, not a cut).",
                 fontsize=10, color="#546775", linespacing=1.6)
        if dropped:
            fig.text(.975, .025, f"Excluded {dropped:,} nonfinite entries", ha="right", fontsize=9)
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        base = output / f"cumulative_scattering_layer_{layer_index:04d}"
        for extension in ["svg", "pdf", "png"]:
            fig.savefig(base.with_suffix("." + extension), dpi=220, facecolor="white")
        print(f"Saved {base}.svg / .pdf / .png")
        print(f"Finite hits: {n}; nonfinite: {dropped}; median: {np.median(angles):.6g} mrad")
        print(f"sigma_core: {sigma:.6g} mrad; full-sample std: {sample_std:.6g} mrad")
        print(f"Beyond +/-3 sigma: Geant4 {tail_fraction:.6g}, Gaussian {gaussian_tail:.6g}")
        if show:
            plt.show()
    return fig, axes


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("root_file")
    parser.add_argument("--layer", type=int, required=True, help="Zero-based layerID")
    parser.add_argument("--output-dir", default="scattering_core")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--sigma-rad", type=float, default=None,
                        help="Optional exact model width, including its low-statistics fallback")
    args = parser.parse_args()
    import uproot
    # Chunked reading avoids loading every layer's angles into memory.
    selected_angles, selected_depths = [], []
    with uproot.open(args.root_file) as root:
        for chunk in root["braggsampler"].iterate(
                ["layerID", "depth", "CumScatteringAngle"], step_size="64 MB", library="np"):
            keep = chunk["layerID"] == args.layer
            if np.any(keep):
                selected_angles.append(chunk["CumScatteringAngle"][keep])
                selected_depths.append(chunk["depth"][keep])
    if not selected_angles:
        raise ValueError(f"No hits found for layer {args.layer}.")
    depths = np.concatenate(selected_depths)
    if not np.all(np.isfinite(depths)) or not np.allclose(depths, depths[0], rtol=1e-7, atol=1e-8):
        raise ValueError("Selected layer has inconsistent downstream depths.")
    plot_scattering_core(np.concatenate(selected_angles), args.layer, depths[0],
                         sigma_rad=args.sigma_rad, output_dir=args.output_dir,
                         show=not args.no_show)


if __name__ == "__main__":
    main()
