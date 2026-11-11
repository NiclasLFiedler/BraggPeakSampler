import sys
import warnings
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import uproot

sys.path.append("../../range_energy/data_analysis")
import analysisFunctions
import mcs_helper as mcs
from scipy.optimize import minimize, minimize_scalar, curve_fit
from scipy.special import expit, erf
from rel_entr_fit_helper import fit_simple_kl, fit_complex_kl
from pathlib import Path
plt.rcParams.update({'font.size': 16})


# Configuration: preserve the settings in the supplied script.
USEPWO = True
USEHIGHLAND = False
USEMODELMIXTURE = True
SAVEDATA = False
SAVEPLOT = True
LARGEOUTPUT = False

ROOT_FILE = "h2oproj.root"
E0 = 220.0
# Set to None to use range_energy(data, E0). These are your current overrides.
R0_OVERRIDE = 6.735 if USEPWO else 30.722  #30.72, 6.73398 6.736 6.728
FRACTION_R0 = 0.8
Z = 73.6 if USEPWO else 7.4 # 73.6, 7.4 due to Mayneord’s formula
X0 = 0.89 if USEPWO else 36.08

N_GEOMETRY = 100_000
N_CONFIGURATIONS = 3000  # Increase and inspect the pointwise sampling error.
CONFIGURATION_SEED = 20120584  # Separate stream from angle/geometry RNG.
# Confirm against your scorer: convert deltaX to cm if it was stored in mm.
DELTA_X_TO_CM = 1.0
PLOT_LAYER_INDEX = 10
X_LIM = 6.0 if USEPWO else 30.3
ANGULAR_MODEL = "mixture"  # "core" or "mixture"
MIXTURE_MODEL = "rutherford" # "rutherford" or "highland"
MIXTURE_FIT_METHOD = "curve_fit"  # "mle", "curve_fit", or "kl"
MIXTURE_BINS = 2000  # Used by curve_fit and the comparison histogram.
CORE_FIT_WIDTH = 2.0

# def gaussian_core_sigma(values, latestSigma):
#     values = values[np.isfinite(values)]
#     if len(values) < 100:
#         return latestSigma*1.01
#     low, high = np.percentile(values, [15.8655, 84.1345])
#     return 0.5 * (high - low)

def gaussian_core_sigma(values, latestSigma=None, fit_width=2.0):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size < 100:
        warnings.warn("Too few angles: reusing the previous core width.")
        return float(latestSigma) if latestSigma is not None else np.nan

    scale = np.diff(np.percentile(values, [15.8655, 84.1345]))[0] / 2
    if scale <= 0:
        return 0.0 if np.all(values == 0) else np.nan

    counts, edges = np.histogram(values / scale, bins=150, range=(-fit_width, fit_width),)
    centers = (edges[:-1] + edges[1:]) / 2

    def gaussian(x, amplitude, sigma):
        return amplitude * np.exp(-0.5 * (x / sigma)**2)

    params, _ = curve_fit(
        gaussian, centers, counts,
        p0=(counts.max(), 1.0),
        sigma=np.sqrt(np.maximum(counts, 1)),
        absolute_sigma=True,
        bounds=([0, 0.05], [np.inf, 5.0]),
    )
    return float(params[1] * scale)

def gaussian_pdf(theta, sigma):
    """Normalized zero-mean Gaussian; theta and sigma use the same units."""
    theta = np.asarray(theta, dtype=float)
    return np.exp(-0.5 * (theta / sigma)**2) / (np.sqrt(2 * np.pi) * sigma)

def gaussian_mixture_pdf(theta, weight, sigma_core, sigma_tail):
    """weight is the integrated probability of the tail component."""
    return (1 - weight) * gaussian_pdf(theta, sigma_core) + (weight) * gaussian_pdf(theta, sigma_tail)

def fit_cumulative_gaussian_mixture(values, min_samples=100, max_iter=5000, tol=1e-6, initial=None):
    """Return (tail_weight, sigma_core, sigma_tail)."""

    angles = np.asarray(values, dtype=float)
    angles = angles[np.isfinite(angles)]

    if angles.size < min_samples:
        raise ValueError("Too few angles to fit.")

    # Normalize for numerical stability; restore angle units at the end.
    scale = np.sqrt(np.mean(angles**2))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Angles must have a finite positive RMS.")

    angle_squared = (angles / scale)**2

    if initial is not None:
        weight, sigma_core, sigma_tail = initial
        weight = 1-weight
        var_core = (sigma_core / scale)**2
        var_tail = (sigma_tail / scale)**2
    else:
        weight = 0.9
        var_core = 0.5**2
        var_tail = 2.0**2

    for _ in range(max_iter):
        log_tail_over_core = ( np.log((1 - weight) / weight) + 0.5 * np.log(var_core / var_tail) + 0.5 * angle_squared * (1 / var_core - 1 / var_tail)
        )
        responsibility = expit(-log_tail_over_core)

        # Update using fractional membership in each component.
        core_count = responsibility.sum()
        tail_count = angles.size - core_count

        if min(core_count, tail_count) <= 0:
            raise RuntimeError("One mixture component became empty.")

        new_weight = core_count / angles.size
        new_var_core = (np.sum(responsibility * angle_squared) / core_count)
        new_var_tail = (np.sum((1 - responsibility) * angle_squared) / tail_count)

        if min(new_var_core, new_var_tail) <= 0:
            raise RuntimeError("A mixture variance collapsed to zero.")

        change = max(
            abs(new_weight - weight),
            abs(np.log(new_var_core / var_core)),
            abs(np.log(new_var_tail / var_tail)),
        )

        weight = new_weight
        var_core = new_var_core
        var_tail = new_var_tail

        if change < tol:
            break
    else:
        raise RuntimeError("Mixture fit did not converge.")

    sigma_core, sigma_tail = scale * np.sqrt([var_core, var_tail])

    if sigma_core > sigma_tail:
        weight = 1 - weight
        sigma_core, sigma_tail = sigma_tail, sigma_core

    return 1-float(weight), float(sigma_core), float(sigma_tail) # return the weight of the tail component, not the core component

def fit_cumulative_gaussian_mixture_curve_fit(values, bins=150, min_samples=100, initial=None):
    angles = np.asarray(values)
    angles = angles[np.isfinite(angles)]

    if len(angles) < min_samples:
        return 0, 0, 0

    counts, edges = np.histogram(angles, bins=bins)
    centres = (edges[:-1] + edges[1:]) / 2
    bin_width = edges[1] - edges[0]

    normalization = len(angles) * bin_width
    density = counts / normalization
    error = np.sqrt(np.maximum(counts, 1)) / normalization

    # Initial guess: 80% narrow Gaussian + 20% broad Gaussian
    sigma = np.std(angles)
    lower_sigma = sigma * 1e-6

    if initial is None:
        initial = [0.05, sigma / 2, 2 * sigma]
    else:
        weight, sigma_core, sigma_tail = initial
        initial = [
            np.clip(weight, 1e-6, 1 - 1e-6),
            max(sigma_core, 2 * lower_sigma),
            max(sigma_tail, 2 * lower_sigma),
        ]

    params, _ = curve_fit(
            gaussian_mixture_pdf,
            centres,
            density,
            p0=initial,
            sigma=error,
            absolute_sigma=True,
            bounds=([0, lower_sigma, lower_sigma],
                    [1, np.inf, np.inf]),
            max_nfev=1000000,
        )

    weight, sigma1, sigma2 = params

    # Always return the narrower component first
    if sigma1 > sigma2:
        weight = 1 - weight
        sigma1, sigma2 = sigma2, sigma1

    return weight, sigma1, sigma2 # returns tail weight

def mixture_mean_nll(values, params):
    """Common in-sample score for either fit; lower is better (same data/units).

    MLE optimizes this score; it is not an independent validation measure.
    Computed outside the fit timer.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    w_tail, s1, s2 = params
    a = np.log1p(-w_tail) - np.log(s1) - 0.5*(x/s1)**2
    b = np.log(w_tail) - np.log(s2) - 0.5*(x/s2)**2
    return float(0.5*np.log(2*np.pi) - np.mean(np.logaddexp(a, b)))


def plot_cumulative_angular_fit(values, core_sigma, layer, depth, mixture=None, alternative=None):
    """Compare the unchanged quantile core approximation with a fitted mixture."""
    angles = np.asarray(values)
    angles = angles[np.isfinite(angles)] * 1000  # rad -> mrad
    if angles.size == 0:
        return
    theta = np.linspace(angles.min(), angles.max(), 2000)
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.hist(angles, bins=MIXTURE_BINS, density=True, histtype="step", linewidth=1.5,
            color="#21618C", label="Geant4")
    ax.plot(theta, gaussian_pdf(theta, core_sigma * 1000), "--", linewidth=2, color="#D55E00", label="Gaussian core (quantile width): {0:.2f} mrad".format(core_sigma * 1000))

    if mixture is not None and np.all(np.isfinite(mixture)):
        w, s1, s2 = mixture
        s1, s2 = s1 * 1000, s2 * 1000
        ax.plot(theta, gaussian_mixture_pdf(theta, w, s1, s2),
                color="black", linewidth=2, label="Mixture: " + ("MLE" if MIXTURE_FIT_METHOD in ("mle", "compare") else MIXTURE_FIT_METHOD))
        ax.plot(theta, (1-w) * gaussian_pdf(theta, s1), ":", color="#009E73",
                label=f"Narrow: w={1-w:.3f}, sigma={s1:.2f} mrad")
        ax.plot(theta, w * gaussian_pdf(theta, s2), ":", color="#CC79A7",
                label=f"Broad: w={w:.3f}, sigma={s2:.2f} mrad")
    if alternative is not None and np.all(np.isfinite(alternative)):
        w, s1, s2 = alternative
        ax.plot(theta, gaussian_mixture_pdf(theta, w, s1*1000, s2*1000),
                "-.", color="#7B3294", linewidth=2, label="Mixture: curve_fit, sigma_1: {0:.2f} mrad, sigma_2: {1:.2f} mrad".format(s1*1000, s2*1000))
    ax.set_yscale("log")
    if np.isfinite(core_sigma) and core_sigma > 0:
        ax.set_ylim(bottom=0.1 / (len(angles) * core_sigma * 1000))
    ax.set_xlabel(r"Cumulative angle $\theta_x$ / mrad")
    ax.set_ylabel(r"Probability density / mrad$^{-1}$")
    ax.set_title(f"Layer {layer} — depth {depth:.2f} cm")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    # fig.savefig(f"scattering_{ANGULAR_MODEL}_layer_{layer}.pdf")
    plt.show()

def integer_ids(values, name):
    values = np.asarray(values)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError(f"{name} must contain finite nonnegative IDs.")
    ids = values.astype(np.int64)
    if not np.all(values == ids):
        raise ValueError(f"{name} contains noninteger values.")
    return ids

def plot_mixture_parameters(z, weight, sigma_core, sigma_tail):
    fig, (ax_weight, ax_sigma) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True
    )

    ax_weight.plot(z, 1 - weight, label="Narrow component", color="#009E73")
    ax_weight.plot(z, weight, label="Broad component", color="#CC79A7")
    ax_weight.set_ylabel("Mixture weight")
    ax_weight.set_ylim(0, 1)
    ax_weight.legend()
    ax_weight.grid(alpha=0.3)

    ax_sigma.plot(
        z, sigma_core * 1000, label="Narrow component", color="#009E73"
    )
    ax_sigma.plot(
        z, sigma_tail * 1000, label="Broad component", color="#CC79A7"
    )
    ax_sigma.set_xlabel("Depth / cm")
    ax_sigma.set_ylabel(r"$\sigma$ / mrad")
    ax_sigma.legend()
    ax_sigma.grid(alpha=0.3)

    fig.tight_layout()
    # fig.savefig("mixture_parameters.pdf")
    plt.show()

def comp_plot_mixture_parameters(z, weight, sigma_core, sigma_tail, z2, weight2, sigma_core2, sigma_tail2):
    fig, (ax_weight, ax_sigma) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True
    )

    ax_weight.plot(z, 1 - weight, label="Narrow component", color="#009E73")
    ax_weight.plot(z, weight, label="Broad component", color="#CC79A7")
    
    ax_weight.plot(z2, 1 - weight2, label="Model: single-kick core probability", color="#FF0000")
    ax_weight.plot(z2, weight2, label="Model: single-kick tail probability", color="#3700FF")
    
    ax_weight.set_ylabel("Mixture weight")
    ax_weight.set_ylim(0, 1)
    ax_weight.legend()
    ax_weight.grid(alpha=0.3)

    ax_sigma.plot(z, sigma_core * 1000, label="Narrow component", color="#009E73")
    ax_sigma.plot(z, sigma_tail * 1000, label="Broad component", color="#CC79A7")
    ax_sigma.plot(z2, sigma_core2 * 1000, label="Model: all-core cumulative width", color="#FF0000")
    ax_sigma.plot(z2, sigma_tail2 * 1000, label="Model: all-tail cumulative width", color="#3700FF")

    ax_sigma.set_xlabel("Depth / cm")
    ax_sigma.set_ylabel(r"$\sigma$ / mrad")
    ax_sigma.legend()
    ax_sigma.grid(alpha=0.3)

    fig.tight_layout()
    # fig.savefig("mixture_parameters.pdf")
    plt.show()

def main():
    if ANGULAR_MODEL not in ("core", "mixture"):
        raise ValueError('ANGULAR_MODEL must be "core" or "mixture".')
    if MIXTURE_FIT_METHOD not in ("mle", "curve_fit", "compare", "kl"):
        raise ValueError('MIXTURE_FIT_METHOD must be "mle", "curve_fit", or "compare".')

    data_file = ("../../range_energy/data_analysis/pbwo4_alt_range_energy.npz"
                 if USEPWO else "../../range_energy/data_analysis/h2o_alt_range_energy.npz")
    data = analysisFunctions.load_EnergyRange(data_file)
    p_exp = float(data.p)
    print(f"Range parameter: {p_exp:.6f}")
    data = analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/pbwo4_range_energy.npz") if USEPWO else analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/h2o_range_energy.npz")
    fitted_R0 = float(analysisFunctions.range_energy(data, E0))
    R0 = fitted_R0 if R0_OVERRIDE is None else float(R0_OVERRIDE)

    print(f"Fitted R0: {fitted_R0:.8f} cm; used R0: {R0:.8f} cm")
    del data, data_file

    with uproot.open(ROOT_FILE) as f:
        tree = f["braggsampler"]
        event = integer_ids(tree["event"].array(library="np"), "event")
        layerID = integer_ids(tree["layerID"].array(library="np"), "layerID")
        # Per-hit coordinates used only to construct z, not a second model grid.
        hit_depth = tree["depth"].array(library="np").astype(float)
        CumScatteringAngle = tree["CumScatteringAngle"].array(library="np")
        ScatteringAngle = tree["SingleScatteringAngle"].array(library="np")
        deltaX = tree["deltaX"].array(library="np") * DELTA_X_TO_CM
        G4ExcessPathLength = tree["ExcessPathLength"].array(library="np")

    if event.size == 0:
        raise ValueError("No recorded hits.")
    if not all(a.shape == event.shape for a in
               (layerID, hit_depth, CumScatteringAngle, ScatteringAngle, deltaX)):
        raise ValueError("ROOT branches must have matching shapes.")

    event_ids, event_index = np.unique(event, return_inverse=True)
    max_layer = np.full(len(event_ids), -1, dtype=np.int64)
    np.maximum.at(max_layer, event_index, layerID)
    last_row = np.r_[event[1:] != event[:-1], True]
    print("Unique recorded events:", len(event_ids))
    print("Contiguous event blocks:", np.count_nonzero(last_row))
    print("Blocks ending below their event's maximum layer:",
          np.count_nonzero(layerID[last_row] < max_layer[event_index[last_row]]))
    n_recorded = len(event_ids)
    del event, event_ids, event_index, last_row

    # Group once by layer ID. This avoids repeatedly scanning all hits and
    # prevents a fixed depth tolerance from mixing neighbouring fine layers.
    order = np.argsort(layerID, kind="stable")
    
    boundary_layers, starts, counts = np.unique(layerID[order], return_index=True, return_counts=True) 

    if not np.array_equal(boundary_layers, np.arange(len(boundary_layers))):
        raise ValueError("Missing layer IDs: need a complete consecutive prefix from 0.")

    # The only depth grid: z[i] = downstream boundary of layer i.
    z = hit_depth[order[starts]].copy()
    if np.any(~np.isfinite(z)) or np.any(z <= 0) or np.any(np.diff(z) <= 0):
        raise ValueError("z must contain finite positive increasing boundaries.")
    n_model = np.searchsorted(z, R0, side="left")
    if n_model == 0:
        raise ValueError("No scored boundary lies below R0.")
    
    CumSigmaCore = np.full(len(z), np.nan)
    SingleSigmaCore = np.full(len(z), np.nan)
    CumDeltaXSigmaCore = np.full(len(z), np.nan)
    # Mixture fits are separate from the original core-width arrays.
    MixtureWeight = np.full(len(z), np.nan)
    CumMixtureSigmaCore = np.full(len(z), np.nan)
    CumMixtureSigmaTail = np.full(len(z), np.nan)

    methods = ("mle", "curve_fit", "kl") if MIXTURE_FIT_METHOD == "compare" else (MIXTURE_FIT_METHOD,)
    fit_parameters = {name: np.full((len(z), 3), np.nan) for name in methods}
    fit_seconds = {name: np.full(len(z), np.nan) for name in methods}
    fit_nll = {name: np.full(len(z), np.nan) for name in methods}
    print(f"Calculating angular widths by layer: {ANGULAR_MODEL}, {MIXTURE_FIT_METHOD}")

    previous_fit = dict.fromkeys(methods)
    previous_width = None
    for i, (start, count) in enumerate(zip(starts[:n_model], counts[:n_model])):
        rows = order[start:start + count]
        recorded_depth = hit_depth[rows]
        if not np.all(np.isfinite(recorded_depth)):
            raise ValueError(f"Nonfinite depths in layer {i}.")
        z[i] = recorded_depth[0]
        if not np.allclose(recorded_depth, z[i], rtol=1e-7, atol=1e-8):
            raise ValueError(f"Layer {i} has inconsistent downstream depths; check scoring.")

        width = z[i] - (z[i-1] if i else 0.0)
        if (previous_width is not None and width < previous_width
                and not np.isclose(width, previous_width, rtol=1e-4, atol=1e-8)):
            previous_fit = dict.fromkeys(methods)
            print(f"Layer {i}: grid spacing {previous_width:.6g} -> {width:.6g} cm; "
                  "resetting mixture fit initial parameters.")
        previous_width = width
        
        CumSigmaCore[i] = gaussian_core_sigma(CumScatteringAngle[rows],
                                              CumSigmaCore[i-1],
                                              #fit_width=CORE_FIT_WIDTH
                                              )
        SingleSigmaCore[i] = gaussian_core_sigma(ScatteringAngle[rows],
                                                 SingleSigmaCore[i-1],
                                                 fit_width=CORE_FIT_WIDTH
                                                 )
        CumDeltaXSigmaCore[i] = gaussian_core_sigma(deltaX[rows],
                                                    CumDeltaXSigmaCore[i-1],
                                                    fit_width=CORE_FIT_WIDTH
                                                    )
        
        mixture = alternative = None

        if ANGULAR_MODEL == "mixture":
            angles = CumScatteringAngle[rows]
            for name in methods:
                start_time = perf_counter()
                try:
                    if name == "mle":
                        params = fit_cumulative_gaussian_mixture(angles, initial=previous_fit[name])
                    elif name == "kl":
                        params = fit_simple_kl(angles, bins=MIXTURE_BINS, initial=previous_fit[name])
                    else:
                        params = fit_cumulative_gaussian_mixture_curve_fit(angles, MIXTURE_BINS, initial=previous_fit[name])
                    previous_fit[name] = params
                    fit_seconds[name][i] = perf_counter() - start_time
                    fit_parameters[name][i] = params
                    fit_nll[name][i] = mixture_mean_nll(angles, params)
                except (ValueError, RuntimeError) as exc:
                    fit_seconds[name][i] = perf_counter() - start_time
                    warnings.warn(f"Layer {i}, {name}: {exc} Parameters remain NaN.")

            selected = "mle" if MIXTURE_FIT_METHOD == "compare" else MIXTURE_FIT_METHOD
            mixture = fit_parameters[selected][i]
            (MixtureWeight[i], CumMixtureSigmaCore[i], CumMixtureSigmaTail[i]) = mixture

            if MIXTURE_FIT_METHOD == "compare":
                alternative = fit_parameters["curve_fit"][i]
        
        if i== PLOT_LAYER_INDEX:
            plot_cumulative_angular_fit(CumScatteringAngle[rows], CumSigmaCore[i], i, z[i], mixture, alternative)
            #plot_cumulative_angular_fit(CumScatteringAngle[rows], CumSigmaCore[i], i, z[i], mixture, alternative)
        if i % 10 == 0:
            print(f"Layer {i}, depth {z[i]:.6f} cm")

    plt.plot(z[:n_model], CumSigmaCore, "o-", markersize=3, label="Cum Core width")

    KickVarCore = SingleSigmaCore[:n_model]**2
    CumVarCore = np.cumsum(KickVarCore)
    CumSigmaCore = np.sqrt(CumVarCore)
    plt.plot(z[:n_model], CumSigmaCore, "-o", markersize=3, label="Cumulative core width from increments ")
    plt.legend()
    plt.show()





    if ANGULAR_MODEL == "mixture":
        CumMixtureVarCore = np.full(len(z), np.nan)
        CumMixtureVarTail = np.full(len(z), np.nan)
        CumMixtureVarCore[:n_model] = mcs.fit_cumulative_variance(CumMixtureSigmaCore[:n_model]**2)
        CumMixtureVarTail[:n_model] = mcs.fit_cumulative_variance(CumMixtureSigmaTail[:n_model]**2)

        CumMixtureSigmaCore = np.sqrt(CumMixtureVarCore)
        CumMixtureSigmaTail = np.sqrt(CumMixtureVarTail)

        MixtureKickVarCore = np.diff(np.r_[0.0, CumMixtureVarCore])
        MixtureKickVarTail = np.diff(np.r_[0.0, CumMixtureVarTail])
        


    if ANGULAR_MODEL == "mixture":
        plt.figure(figsize=(10, 7))
        plt.plot(z, CumMixtureVarCore, "o-", markersize=3,
                 label="Core variance")
        plt.plot(z, CumMixtureVarTail, "^-", markersize=3,
                 label="Tail variance")

        plt.xlabel("Depth / cm")
        plt.ylabel(r"Cumulative angular variance / rad$^2$")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.show()
        CumMixtureVariance = ((1-MixtureWeight) * CumMixtureVarCore + (MixtureWeight) * CumMixtureVarTail)
        if SAVEDATA:
            np.savez("cumulative_angular_mixture.npz", z=z,
                 weight_core=1-MixtureWeight, weight_tail=MixtureWeight,
                 sigma_core=np.sqrt(CumMixtureVarCore), sigma_tail=np.sqrt(CumMixtureVarTail),
                 variance=CumMixtureVariance)
            print("Saved mixture fits (angles in rad): cumulative_angular_mixture.npz")
            
    del order, starts, counts, hit_depth
    del CumScatteringAngle, ScatteringAngle, deltaX
    
    if np.any(z <= 0) or np.any(np.diff(z) <= 0):
        raise ValueError("z must contain positive increasing downstream boundaries.")
    # Include the interval from the material entrance (depth 0) to z[0].
    dz = np.diff(np.r_[0.0, z])
    print(f"Intervals: {len(z)}, dz min/max: {dz.min():.6g}/{dz.max():.6g} cm")

    max_layer.sort()
    n_reaching = len(max_layer) - np.searchsorted(max_layer, boundary_layers, side="left")
    n_primaries = n_recorded
    if not np.isfinite(n_primaries) or n_primaries < n_recorded or int(n_primaries) != n_primaries:
        raise ValueError("N_PRIMARIES must be an integer >= recorded event count.")
    reach_g4 = n_reaching / n_primaries
    
    del max_layer
    
    CumVarCore = CumSigmaCore**2  # Input angles are already radians.
    CumVarHighland = None
    CumMixtureModelWeight = np.full(len(z), np.nan)
    CumMixtureModelSigmaCore = np.full(len(z), np.nan)
    CumMixtureModelSigmaTail = np.full(len(z), np.nan)
    CumMixtureModelVarCore = np.full(len(z), np.nan)
    CumMixtureModelVarTail = np.full(len(z), np.nan)
    
    
    CumVarCoreFit = mcs.fit_cumulative_variance(CumVarCore[:n_model])
    V = CumVarCoreFit
    CumVarHighland = mcs.highland_variance(z=z[:n_model], dz=dz[:n_model], R0=R0, E0=E0, p_exp=p_exp, X0=X0)
    if USEHIGHLAND:
        V = CumVarHighland
    print(V)
    KickVarCore = np.diff(np.r_[0.0, V])
    
    Xs = mcs.scattering_length_from_X0(X0, Z)
    MixtureModelKickVarCore, MixtureModelKickVarTail, weight = mcs.gaussian_mixture_model(z=z[:n_model],
                                                                                          dz=dz[:n_model],
                                                                                          R0=R0,
                                                                                          E0=E0,
                                                                                          p_exp=p_exp,
                                                                                          X0=X0,
                                                                                          Z=Z,
                                                                                          Xs=Xs,
                                                                                          KickVariance=KickVarCore[:n_model],
                                                                                          MIXTURE_MODEL=MIXTURE_MODEL,)
    CumMixtureModelWeight = weight
    CumMixtureModelVarCore = np.cumsum(MixtureModelKickVarCore)
    CumMixtureModelVarTail = np.cumsum(MixtureModelKickVarTail)
    CumMixtureModelSigmaCore = np.sqrt(CumMixtureModelVarCore)
    CumMixtureModelSigmaTail = np.sqrt(CumMixtureModelVarTail)

    if USEMODELMIXTURE:
        comp_plot_mixture_parameters(
            z=z,
            weight=MixtureWeight,
            sigma_core=CumMixtureSigmaCore,
            sigma_tail=CumMixtureSigmaTail,
            z2=z[:n_model],
            weight2=CumMixtureModelWeight[:n_model],
            sigma_core2=CumMixtureModelSigmaCore[:n_model],
            sigma_tail2=CumMixtureModelSigmaTail[:n_model]
        )
        MixtureWeight = CumMixtureModelWeight
        CumMixtureSigmaCore = CumMixtureModelSigmaCore
        CumMixtureSigmaTail = CumMixtureModelSigmaTail
        MixtureKickVarCore = MixtureModelKickVarCore
        MixtureKickVarTail = MixtureModelKickVarTail
    
    plt.rcParams.update({
        "font.size": 16,
        "axes.labelsize": 11,
        "legend.fontsize": 8,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9    ,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.prop_cycle": plt.cycler(
            color=["#000000", "#0072B2", "#D55E00"]
        ),
    })

    g4color = "#000000"
    singularColor ="#0072B2" 
    mixtureColor = "#D55E00"

    fig, (ax, ax_res) = plt.subplots(
        2, 1,
        sharex=True,
        figsize=(15 / 2.54, 14 / 2.54),
        gridspec_kw={"height_ratios": [3, 1]},
        constrained_layout=True,
    )

    target_depth = FRACTION_R0 * R0
    layer_index = int(np.argmin(np.abs(z[:n_model] - target_depth)))

    options = dict(
        layers=[layer_index],
        x_max_mm=2.0,
        bins=150,
        n_geometry=N_GEOMETRY,
        R0=R0,
        geant4=(layerID, G4ExcessPathLength),
        ax=ax,
        ax_res=ax_res,
    )
    print("Core")
    mcs.plot_excess_pathlength_layers(
        z, dz,
        cumulative_variance=V,
        **options,
        Range=R0
    )
    del layerID, G4ExcessPathLength

    if ANGULAR_MODEL== "mixture":
        print("Anuglar Model")
        mcs.plot_excess_pathlength_layers(
            z, dz,
            mixture=(MixtureKickVarCore[:n_model], MixtureKickVarTail[:n_model], MixtureWeight[:n_model]),
            show_geant4=False,
            **options,
            Range=R0
        )
    
    ax.set_xlabel("")
    ax_res.set(
        xlabel=r"Excess path length $\Delta r$ / mm",
        ylabel="Residual\n" + r"/ mm$^{-1}$",
    )
    ax.set_axisbelow(True)
    ax.grid(
        True, which="major",
        color="0.9", linestyle="-", linewidth=0.7,
    )
    ax_res.axhline(0, color="0.4", linewidth=0.8, zorder=0)
    ax_res.grid(
        True, which="major",
        color="0.9", linestyle="-", linewidth=0.7,
    )
    if SAVEPLOT:
        if USEPWO:
            fig.savefig(f"img/pwo/excess_pathlength_layers_{E0}.pdf", format="pdf", dpi=300)
        else:
            fig.savefig(f"img/h2o/excess_pathlength_layers_{E0}.pdf", format="pdf", dpi=300)
    plt.show()
    del options

    z_fit = z[:n_model]
    plt.figure(figsize=(12, 9))

    if CumVarHighland is not None:
        plt.plot(z_fit, np.degrees(np.sqrt(CumVarHighland[:n_model])), label="Integrated Highland, global log",)

    plt.scatter(z_fit, np.degrees(CumSigmaCore[:n_model]), s=10, label="Geant4 cumulative core width",)

    if ANGULAR_MODEL == "mixture":
        plt.plot(z_fit, np.degrees(CumMixtureSigmaCore[:n_model]), label="Mixture narrow sigma",)
        plt.plot(z_fit, np.degrees(CumMixtureSigmaTail[:n_model]), label="Mixture broad sigma",)
        plt.plot(z_fit, np.degrees(np.sqrt(CumMixtureVariance[:n_model])), "--", label="Mixture RMS",)

        
        var_core = CumMixtureModelVarCore[:n_model]
        var_tail = CumMixtureModelVarTail[:n_model]
        weight_tail = CumMixtureModelWeight[:n_model]
        var_total = (1 - weight_tail) * var_core + weight_tail * var_tail

        plt.plot(z_fit, np.degrees(np.sqrt(var_core)),":", label="Mixture model narrow sigma",)
        plt.plot(z_fit, np.degrees(np.sqrt(var_tail)), ":", label="Mixture model broad sigma",)
        plt.plot(z_fit, np.degrees(np.sqrt(var_total)), "--", color="orange",label="Mixture model RMS",)

    plt.xlabel("Depth / cm")
    plt.ylabel("Projected angle / degree")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

    if USEHIGHLAND:
        mcs.simulation_grid(z[:n_model], V[:n_model], R0, layer_ids=boundary_layers[:n_model])

    rng = np.random.default_rng(12345)
    reach_analytical = np.zeros(len(z))

    reach_geometry = np.zeros(len(z))
    reach_geometry_error = np.zeros(len(z))

    mixture_reach_analytical = np.zeros(len(z))
    mixture_reach_analytical_error = np.zeros(len(z))
    mixture_reach_geometry = np.zeros(len(z))
    mixture_reach_geometry_error = np.zeros(len(z))
    
    if ANGULAR_MODEL == "mixture":
        def report_configuration_progress(done, total, unique):
            if done == 1 or done % 10 == 0 or done == total:
                print(f"Sampled mixture: depth {done}/{total}, "
                      f"{unique} unique prefixes from {N_CONFIGURATIONS} configurations")
        print(f"Length: {z[:n_model]}, {dz[:n_model]}, R0: {R0}, "
              f"KickVarCore: {MixtureKickVarCore[:n_model]}, "
              f"KickVarTail: {MixtureKickVarTail[:n_model]}, "
              f"Weight: {MixtureWeight[:n_model]}")
        (mixture_reach_analytical[:n_model], mixture_reach_analytical_error[:n_model]) = mcs.reach_probability_sampled_mixture(z[:n_model],
                                                                                                           dz[:n_model],
                                                                                                           R0,
                                                                                                           MixtureKickVarCore[:n_model],
                                                                                                           MixtureKickVarTail[:n_model],
                                                                                                           MixtureWeight[:n_model],
                                                                                                           np.random.default_rng(CONFIGURATION_SEED),
                                                                                                           n_configurations=N_CONFIGURATIONS,
                                                                                                           progress=report_configuration_progress,
                                                                                                           )
        for j in range(n_model):
            (mixture_reach_geometry[j], mixture_reach_geometry_error[j]) = mcs.reach_probability_rng_geometry_mixture(R0-z[j],
                                                                                                                     MixtureWeight[:j+1],
                                                                                                                     MixtureKickVarCore[:j+1],
                                                                                                                     MixtureKickVarTail[:j+1],
                                                                                                                     dz[:j+1],
                                                                                                                     rng,
                                                                                                                     n_samples=N_GEOMETRY,
                                                                                                                     )
        print("Maximum configuration-sampling standard error: " 
              f"{mixture_reach_analytical_error[:n_model].max():.4g}")
        if SAVEDATA:
            np.savez("sampled_mixture_reach.npz", 
                     z=z[:n_model],
                     reach=mixture_reach_analytical[:n_model],
                     standard_error=mixture_reach_analytical_error[:n_model],
                     n_configurations=N_CONFIGURATIONS,
                     seed=CONFIGURATION_SEED,
                     var_core=MixtureKickVarCore[:n_model],
                     var_tail=MixtureKickVarTail[:n_model],
                     p_tail=MixtureWeight[:n_model],
                     dz=dz[:n_model],
                     R0=R0,
                     )

    for j in range(n_model):
        remaining_range = R0 - z[j]
        V_j, dz_j = V[:j+1], dz[:j+1]
        weights = mcs.weighted_eigenvalues(V_j, dz_j)
        probability = float(mcs.gchi2_exact_cdf(remaining_range, weights))
        reach_analytical[j] = probability
        reach_geometry[j], reach_geometry_error[j] = mcs.reach_probability_rng_geometry(remaining_range, V_j, dz_j, rng, n_samples=N_GEOMETRY)

    z_plot = np.r_[z, z[-1] + dz[-1]]
    dz_plot = np.r_[dz, dz[-1]]

    reach_g4_plot = np.r_[reach_g4, 0.0]
    reach_analytical_plot = np.r_[reach_analytical, 0.0]
    reach_geometry_plot = np.r_[reach_geometry, 0.0]

    reach_geometry_error_plot = np.r_[reach_geometry_error, 0.0]

    mixture_reach_analytical_plot = np.r_[mixture_reach_analytical, 0.0]
    mixture_reach_geometry_plot = np.r_[mixture_reach_geometry, 0.0]

    mixture_reach_geometry_error_plot = np.r_[mixture_reach_geometry_error, 0.0]
    
    fig_reach, fig_stop = mcs.plot_reach_and_stopping(
        z_plot,
        reach_g4_plot,
        reach_analytical_plot,
        reach_geometry_plot=reach_geometry_plot,
        mixture_reach_analytical_plot=(
            mixture_reach_analytical_plot if ANGULAR_MODEL== "mixture" else None
        ),
        mixture_reach_geometry_plot=(
            mixture_reach_geometry_plot if ANGULAR_MODEL== "mixture" else None
        ),
        R0=R0,
        n_primaries=n_primaries,
        n_configurations=N_CONFIGURATIONS,
        xlim=(X_LIM, z_plot[-1] + 0.03),
    )

    if SAVEPLOT:
        output_dir = Path("img") / ("pwo" if USEPWO else "h2o")
        output_dir.mkdir(parents=True, exist_ok=True)

        with plt.rc_context({"pdf.fonttype": 42}):
            fig_reach.savefig(output_dir / f"reach_probability_{E0}.pdf")
            fig_stop.savefig(output_dir / f"stopping_density_{E0}.pdf")

    plt.show()


    print(f"Last model evaluation: {z[n_model-1]:.8f} cm; "
          f"gap to R0: {R0-z[n_model-1]:.8f} cm; "
          f"reach there: {reach_analytical[n_model-1]:.6g}")
    print(f"Geant4 reach at last recorded boundary: {reach_g4[-1]:.6g}")

    # Diagnostic increments: leave negative estimates undefined rather than
    # silently clipping them to zero. They can arise from noisy/core widths.
    increment_variance = np.diff(np.r_[0.0, CumVarCore])
    SingleAngleRMSFromCum = np.full(len(z), np.nan)
    good = np.isfinite(increment_variance) & (increment_variance >= 0)
    SingleAngleRMSFromCum[good] = np.sqrt(increment_variance[good])

    plt.figure(figsize=(12, 9))
    plt.plot(z, np.degrees(SingleAngleRMSFromCum), "o-", label="From cumulative core variance")
    plt.plot(z, np.degrees(SingleSigmaCore), "s--", label="Direct Geant4 single-angle core width")
    plt.xlabel("Depth / cm")
    plt.ylabel("Single projected angle / degree")
    plt.legend()
    plt.grid(alpha=.3)
    plt.tight_layout()
    plt.show()

    # Retain your lever-arm diagnostic on unequal intervals. With right-end 
    # angles held over each interval, kick i acts over z[j]-(z[i]-dz[i]).
    # This is a SMALL-ANGLE variance estimate: use Var(theta), not tan(RMS)^2.
    lateralVariance = np.full(len(z), np.nan)
    for j, d in enumerate(z):
        if np.all(good[:j+1]):
            lever_arm = d-z[:j+1]+dz[:j+1]
            lateralVariance[j] = np.sum(lever_arm**2 * increment_variance[:j+1])
    plt.figure(figsize=(12, 9))
    plt.plot(z, np.sqrt(lateralVariance), "o-", label="Gaussian increment model")
    plt.plot(z, CumDeltaXSigmaCore, label="Geant4 lateral core width")
    plt.xlabel("Depth / cm")
    plt.ylabel("Lateral width / cm")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    # plt.savefig("lateral_scattering.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()