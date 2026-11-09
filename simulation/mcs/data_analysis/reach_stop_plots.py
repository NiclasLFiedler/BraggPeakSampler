import matplotlib.pyplot as plt
import numpy as np


def plot_reach_and_stopping(
    z_plot, reach_g4_plot, reach_analytical_plot, *,
    reach_geometry_plot=None, mixture_reach_analytical_plot=None,
    mixture_reach_geometry_plot=None, R0=None, n_primaries=None, xlim=None,
):
    """Return reach and stopping figures; depths in cm, residuals model - Geant4.
    Inputs share an increasing depth grid. Geant4 reach must be nonincreasing.
    Optional errors assume n_primaries independent, equally weighted primaries.
    """
    z = np.asarray(z_plot)
    g4 = np.asarray(reach_g4_plot)
    widths = np.diff(z)
    models = [
        ("Single Gaussian", reach_analytical_plot, "#0072B2", "-"),
        ("Single Gaussian, sampled", reach_geometry_plot, "#0072B2", "--"),
        ("Gaussian mixture", mixture_reach_analytical_plot, "#D55E00", "-."),
        ("Gaussian mixture, sampled", mixture_reach_geometry_plot, "#D55E00", ":"),
    ]

    figures = []
    with plt.rc_context({"font.size": 9, "legend.fontsize": 8}):
        for stopping in (False, True):
            fig, (ax, ax_res) = plt.subplots(
                2, 1, sharex=True, figsize=(15 / 2.54, 14 / 2.54),
                gridspec_kw={"height_ratios": [3, 1]}, constrained_layout=True,
            )
            x = 0.5 * (z[:-1] + z[1:]) if stopping else z
            probability = -np.diff(g4) if stopping else g4
            reference = probability / widths if stopping else probability
            error = None
            if n_primaries is not None:
                error = np.sqrt(probability * (1 - probability) / n_primaries)
                if stopping:
                    error = error / widths
                band = dict(color="0.85", linewidth=0, zorder=0,
                            label=r"Geant4 stat. uncertainty ($\pm1\sigma$)")
                if stopping:
                    ax_res.stairs(error, z, baseline=-error, fill=True, **band)
                else:
                    ax_res.fill_between(z, -error, error, **band)
                
            ax.errorbar(
                x, reference, yerr=error, fmt="o", markersize=3,
                markerfacecolor="none", color="black", elinewidth=0.8,
                capsize=0, zorder=3, label="Geant4",
            )
            for label, reach, color, linestyle in models:
                if reach is None:
                    continue
                line = dict(color=color, linestyle=linestyle, linewidth=1.4)
                if stopping:
                    density = -np.diff(reach) / widths
                    ax.stairs(density, z, baseline=None, label=label, **line)
                    ax_res.stairs(density - reference, z, baseline=None, **line)
                else:
                    ax.plot(z, reach, label=label, **line)
                    ax_res.plot(z, np.asarray(reach) - g4, **line)

            for panel in (ax, ax_res):
                panel.grid(color="0.9", linewidth=0.5)
                if R0 is not None:
                    panel.axvline(R0, color="0.5", linestyle="--", linewidth=1,
                                  label=r"$ \mathrm{R_0}$" if panel is ax else "_nolegend_")
            ax_res.axhline(0, color="0.4", linewidth=0.8)
            ax_res.set_xlabel(r"Depth $z$ / cm")
            ax_res.set_ylabel(r"Model $-$ Geant4" +
                              ("\n" + r"/ cm$^{-1}$" if stopping else ""))
            ax.set_ylabel(r"Stopping density / cm$^{-1}$" if stopping
                          else "Reach probability")
            if not stopping:
                ax.set_ylim(-0.03, 1.03)
            ax_res.set_xlim(xlim if xlim is not None else (z[0], z[-1]))
            ax.legend()
            figures.append(fig)
    return tuple(figures)
