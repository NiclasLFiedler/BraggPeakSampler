import sys
import warnings
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import uproot

sys.path.append("../../range_energy/data_analysis")
import analysisFunctions
import mcs_helper as mcs
from plotter import plot_scattering_core
from scipy.optimize import minimize, curve_fit
from scipy.special import expit

plt.rcParams.update({'font.size': 16})

# Configuration: preserve the settings in the supplied script.
USEPWO = False
USEHIGHLAND = False
USEMODELMIXTURE = True
SAVEDATA = False
LARGEOUTPUT = False

ROOT_FILE = "h2oproj.root"
E0 = 220.0
# Set to None to use range_energy(data, E0). These are your current overrides.
R0_OVERRIDE = 6.74 if USEPWO else 30.75  #30.72, 6.73398 6.736 6.728
Z = 73.6 if USEPWO else 7.4 # 73.6, 7.4 due to Mayneord’s formula
X0 = 0.89 if USEPWO else 36.08
# Set to the generated-primary count to include events with no recorded hit.
# None preserves your original normalization to recorded events only.
N_PRIMARIES = None
N_RNG = 100_000
N_GEOMETRY = 100_000
N_CONFIGURATIONS = 10000  # Increase and inspect the pointwise sampling error.
CONFIGURATION_SEED = 20120584  # Separate stream from angle/geometry RNG.
# Confirm against your scorer: convert deltaX to cm if it was stored in mm.
DELTA_X_TO_CM = 1.0
PLOT_LAYER_INDEX = 60
X_LIM = 6.0 if USEPWO else 30.3
ANGULAR_MODEL = "core"  # "core" or "mixture"
MIXTURE_MODEL = "rutherford" # "rutherford" or "highland"
MIXTURE_FIT_METHOD = "curve_fit"  # "mle", "curve_fit", or "compare"
MIXTURE_BINS = 500  # Used by curve_fit and the comparison histogram.

def gaussian_core_sigma(values, latestSigma):
    values = values[np.isfinite(values)]
    if len(values) < 100:
        return latestSigma*1.01
    low, high = np.percentile(values, [15.8655, 84.1345])
    return 0.5 * (high - low)

def gaussian_pdf(theta, sigma):
    """Normalized zero-mean Gaussian; theta and sigma use the same units."""
    theta = np.asarray(theta, dtype=float)
    return np.exp(-0.5 * (theta / sigma)**2) / (np.sqrt(2 * np.pi) * sigma)


def gaussian_mixture_pdf(theta, weight, sigma_core, sigma_tail):
    """weight is the integrated probability of the narrower component."""
    return (weight * gaussian_pdf(theta, sigma_core)
            + (1 - weight) * gaussian_pdf(theta, sigma_tail))

def fit_cumulative_gaussian_mixture(values, min_samples=100, max_iter=5000, tol=1e-6, initial=None):
    """Return (core_weight, sigma_core, sigma_tail)."""

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

    return float(weight), float(sigma_core), float(sigma_tail)

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
    initial = [0.90, sigma / 2, 2 * sigma]

    lower_sigma = sigma * 1e-6

    if initial is None:
        initial = [0.9, sigma / 2, 2 * sigma]
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

    return weight, sigma1, sigma2

def mixture_mean_nll(values, params):
    """Common in-sample score for either fit; lower is better (same data/units).

    MLE optimizes this score; it is not an independent validation measure.
    Computed outside the fit timer.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    w, s1, s2 = params
    a = np.log(w) - np.log(s1) - 0.5*(x/s1)**2
    b = np.log1p(-w) - np.log(s2) - 0.5*(x/s2)**2
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
                color="black", linewidth=2, label="Mixture: " + ("curve_fit" if MIXTURE_FIT_METHOD == "curve_fit" else "MLE"))
        ax.plot(theta, w * gaussian_pdf(theta, s1), ":", color="#009E73",
                label=f"Narrow: w={w:.3f}, sigma={s1:.2f} mrad")
        ax.plot(theta, (1-w) * gaussian_pdf(theta, s2), ":", color="#CC79A7",
                label=f"Broad: w={1-w:.3f}, sigma={s2:.2f} mrad")
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
    fig.savefig(f"scattering_{ANGULAR_MODEL}_layer_{layer}.svg")
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

    ax_weight.plot(z, weight, label="Narrow component", color="#009E73")
    ax_weight.plot(z, 1 - weight, label="Broad component", color="#CC79A7")
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
    fig.savefig("mixture_parameters.svg")
    plt.show()

def comp_plot_mixture_parameters(z, weight, sigma_core, sigma_tail, z2, weight2, sigma_core2, sigma_tail2):
    fig, (ax_weight, ax_sigma) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True
    )

    ax_weight.plot(z, weight, label="Narrow component", color="#009E73")
    ax_weight.plot(z, 1 - weight, label="Broad component", color="#CC79A7")
    
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
    fig.savefig("mixture_parameters.svg")
    plt.show()

def main():
    if ANGULAR_MODEL not in ("core", "mixture"):
        raise ValueError('ANGULAR_MODEL must be "core" or "mixture".')
    if MIXTURE_FIT_METHOD not in ("mle", "curve_fit", "compare"):
        raise ValueError('MIXTURE_FIT_METHOD must be "mle", "curve_fit", or "compare".')
    with uproot.open(ROOT_FILE) as f:
        tree = f["braggsampler"]
        event = integer_ids(tree["event"].array(library="np"), "event")
        layerID = integer_ids(tree["layerID"].array(library="np"), "layerID")
        # Per-hit coordinates used only to construct z, not a second model grid.
        hit_depth = tree["depth"].array(library="np").astype(float)
        CumScatteringAngle = tree["CumScatteringAngle"].array(library="np")
        ScatteringAngle = tree["SingleScatteringAngle"].array(library="np")
        deltaX = tree["deltaX"].array(library="np") * DELTA_X_TO_CM

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
    z = np.empty(len(boundary_layers), dtype=float)
    CumSigmaCore = np.full(len(z), np.nan)
    SingleSigmaCore = np.full(len(z), np.nan)
    CumDeltaXSigmaCore = np.full(len(z), np.nan)
    # Mixture fits are separate from the original core-width arrays.
    MixtureWeight = np.full(len(z), np.nan)
    CumMixtureSigmaCore = np.full(len(z), np.nan)
    CumMixtureSigmaTail = np.full(len(z), np.nan)

    methods = ("mle", "curve_fit") if MIXTURE_FIT_METHOD == "compare" else (MIXTURE_FIT_METHOD,)
    fit_parameters = {name: np.full((len(z), 3), np.nan) for name in methods}
    fit_seconds = {name: np.full(len(z), np.nan) for name in methods}
    fit_nll = {name: np.full(len(z), np.nan) for name in methods}
    print(f"Calculating angular widths by layer: {ANGULAR_MODEL}, {MIXTURE_FIT_METHOD}")

    previous_curve_fit = None
    for i, (start, count) in enumerate(zip(starts, counts)):
        rows = order[start:start + count]
        recorded_depth = hit_depth[rows]
        if not np.all(np.isfinite(recorded_depth)):
            raise ValueError(f"Nonfinite depths in layer {i}.")
        z[i] = recorded_depth[0]
        if not np.allclose(recorded_depth, z[i], rtol=1e-7, atol=1e-8):
            raise ValueError(f"Layer {i} has inconsistent downstream depths; check scoring.")
        CumSigmaCore[i] = gaussian_core_sigma(CumScatteringAngle[rows], CumSigmaCore[i-1])
        SingleSigmaCore[i] = gaussian_core_sigma(ScatteringAngle[rows], SingleSigmaCore[i-1])
        CumDeltaXSigmaCore[i] = gaussian_core_sigma(deltaX[rows], CumDeltaXSigmaCore[i-1])
        
        mixture = alternative = None

        if ANGULAR_MODEL == "mixture":
            angles = CumScatteringAngle[rows]
            for name in methods:
                start_time = perf_counter()
                try:
                    if name == "mle":
                        params = fit_cumulative_gaussian_mixture(angles, initial=previous_curve_fit)
                        previous_curve_fit = params
                    else:
                        params = fit_cumulative_gaussian_mixture_curve_fit(angles, MIXTURE_BINS, initial=previous_curve_fit,)
                        previous_curve_fit = params
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

        MixtureKickVarCore = np.diff(np.r_[0.0, CumMixtureSigmaCore**2])
        MixtureKickVarTail = np.diff(np.r_[0.0, CumMixtureSigmaTail**2])
        
        if i== PLOT_LAYER_INDEX:
            plot_cumulative_angular_fit(CumScatteringAngle[rows], CumSigmaCore[i], i, z[i], mixture, alternative)
            #plot_cumulative_angular_fit(CumScatteringAngle[rows], CumSigmaCore[i], i, z[i], mixture, alternative)
        if i % 10 == 0:
            print(f"Layer {i}, depth {z[i]:.6f} cm")

    plt.figure(figsize=(10, 7))
    plt.plot(z, CumMixtureSigmaCore**2, "o-", markersize=3,
             label="Core variance")
    plt.plot(z, CumMixtureSigmaTail**2, "^-", markersize=3,
             label="Tail variance")

    plt.xlabel("Depth / cm")
    plt.ylabel(r"Cumulative angular variance / rad$^2$")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
    if ANGULAR_MODEL == "mixture":
        # Second moment of the fitted marginal distribution, not a Gaussian core.
        CumMixtureVariance = (MixtureWeight * CumMixtureSigmaCore**2
                              + (1-MixtureWeight) * CumMixtureSigmaTail**2)
        if SAVEDATA:
            np.savez("cumulative_angular_mixture.npz", z=z,
                 weight_core=MixtureWeight, weight_tail=1-MixtureWeight,
                 sigma_core=CumMixtureSigmaCore, sigma_tail=CumMixtureSigmaTail,
                 variance=CumMixtureVariance)
            print("Saved mixture fits (angles in rad): cumulative_angular_mixture.npz")
            print("Reach uses sampled single-kick mixtures when USEHIGHLAND is True. "
                  "Per-depth cumulative mixture fits are diagnostics only.")
                
    del order, starts, counts, layerID, hit_depth
    del CumScatteringAngle, ScatteringAngle, deltaX
    
    if np.any(z <= 0) or np.any(np.diff(z) <= 0):
        raise ValueError("z must contain positive increasing downstream boundaries.")
    # Include the interval from the material entrance (depth 0) to z[0].
    dz = np.diff(np.r_[0.0, z])
    print(f"Intervals: {len(z)}, dz min/max: {dz.min():.6g}/{dz.max():.6g} cm")

    # Geant4 reach at exactly the same z entries. Do not invent an extra
    # boundary or infer its width from the preceding layer.
    max_layer.sort()
    n_reaching = len(max_layer) - np.searchsorted(max_layer, boundary_layers, side="left")
    n_primaries = n_recorded if N_PRIMARIES is None else N_PRIMARIES
    if not np.isfinite(n_primaries) or n_primaries < n_recorded or int(n_primaries) != n_primaries:
        raise ValueError("N_PRIMARIES must be an integer >= recorded event count.")
    reach_analytical_g4 = n_reaching / n_primaries
    
    del max_layer
    
    if N_PRIMARIES is None:
        print("Geant4 reach is normalized to recorded events (N_PRIMARIES=None).")

    data_file = ("../../range_energy/data_analysis/pbwo4_alt_range_energy.npz"
                 if USEPWO else "../../range_energy/data_analysis/h2o_alt_range_energy.npz")
    data = analysisFunctions.load_EnergyRange(data_file)
    p_exp = float(data.p)
    data = analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/pbwo4_range_energy.npz") if USEPWO else analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/h2o_range_energy.npz")
    fitted_R0 = float(analysisFunctions.range_energy(data, E0))
    R0 = fitted_R0 if R0_OVERRIDE is None else float(R0_OVERRIDE)
    
    if not np.isfinite(R0) or R0 <= 0:
        raise ValueError("R0 must be finite and positive.")
    print(f"Fitted R0: {fitted_R0:.8f} cm; used R0: {R0:.8f} cm")

    # Keep z and every per-depth array at their full original length.
    # The fixed-range model is evaluated only on this contiguous prefix.
    
    n_model = np.searchsorted(z, R0, side="left") #use left
    if n_model == 0:
        raise ValueError("No scored boundary lies below R0.")
    
    CumVarCore = CumSigmaCore**2  # Input angles are already radians.
    CumVarHighland = np.full(len(z), np.nan)
    CumMixtureModelWeight = np.full(len(z), np.nan)
    CumMixtureModelSigmaCore = np.full(len(z), np.nan)
    CumMixtureModelSigmaTail = np.full(len(z), np.nan)
    CumMixtureModelVarCore = np.full(len(z), np.nan)
    CumMixtureModelVarTail = np.full(len(z), np.nan)
    
    CumVarHighland[:n_model] = mcs.highland_variance(z=z[:n_model], dz=dz[:n_model], R0=R0, E0=E0, p_exp=p_exp, X0=X0)
    
    KickVarCore = np.diff(np.r_[0.0, CumVarCore]) if not USEHIGHLAND else np.diff(np.r_[0.0, CumVarHighland])
    KickVarCore = np.diff(np.r_[0.0, CumVarCore])

    if ANGULAR_MODEL == "mixture":
        Xs = mcs.scattering_length_from_X0(X0, Z)
        Xs = None
        var1, var2, weight = mcs.gaussian_mixture_model(z=z[:n_model], dz=dz[:n_model], R0=R0, E0=E0, p_exp=p_exp, X0=X0, Z=Z, Xs=Xs, KickVariance=KickVarCore[:n_model], MIXTURE_MODEL=MIXTURE_MODEL)
        CumMixtureModelWeight = weight
        MixtureModelKickVarCore = var1
        MixtureModelKickVarTail = var2
        CumMixtureModelVarCore = np.cumsum(var1)
        CumMixtureModelVarTail = np.cumsum(var2)
        CumMixtureModelSigmaCore = np.sqrt(CumMixtureModelVarCore)
        CumMixtureModelSigmaTail = np.sqrt(CumMixtureModelVarTail)

    if USEMODELMIXTURE and ANGULAR_MODEL == "mixture":
        comp_plot_mixture_parameters(
            z=z,
            weight=MixtureWeight,
            sigma_core=CumMixtureSigmaCore,
            sigma_tail=CumMixtureSigmaTail,
            z2=z[:n_model],
            weight2=CumMixtureModelWeight,
            sigma_core2=CumMixtureModelSigmaCore,
            sigma_tail2=CumMixtureModelSigmaTail
        )
        MixtureWeight = CumMixtureModelWeight
        CumMixtureSigmaCore = CumMixtureModelSigmaCore
        CumMixtureSigmaTail = CumMixtureModelSigmaTail
        MixtureKickVarCore = MixtureModelKickVarCore
        MixtureKickVarTail = MixtureModelKickVarTail

    
    V = CumVarHighland if USEHIGHLAND else CumVarCore
    mcs.simulation_grid(z[:n_model], V[:n_model], R0,
                        layer_ids=boundary_layers[:n_model])

    rng = np.random.default_rng(12345)
    reach_analytical = np.zeros(len(z))
    stop_analytical = np.zeros(len(z))
    reach_rng = np.zeros(len(z))
    reach_rng_error = np.zeros(len(z))
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
        (mixture_reach_analytical, mixture_reach_analytical_error) = mcs.reach_probability_sampled_mixture(
            z[:n_model], dz[:n_model], R0, MixtureKickVarCore[:n_model], MixtureKickVarTail[:n_model], MixtureWeight[:n_model],
            np.random.default_rng(CONFIGURATION_SEED),
            n_configurations=N_CONFIGURATIONS,
            progress=report_configuration_progress)
        print("Maximum configuration-sampling standard error: "
              f"{mixture_reach_analytical_error[:n_model].max():.4g}")
        if SAVEDATA:
            np.savez("sampled_mixture_reach.npz", z=z[:n_model],
                     reach=mixture_reach_analytical[:n_model],
                     standard_error=mixture_reach_analytical_error[:n_model],
                     n_configurations=N_CONFIGURATIONS, seed=CONFIGURATION_SEED,
                     var_core=MixtureKickVarCore[:n_model], var_tail=MixtureKickVarTail[:n_model], p_tail=MixtureWeight[:n_model],
                     dz=dz[:n_model], R0=R0)
        for j in range(n_model):
            mixture_reach_geometry[j], mixture_reach_geometry_error[j] = (
                mcs.reach_probability_rng_geometry_mixture(
                    R0-z[j], MixtureWeight[:j+1],
                    MixtureKickVarCore[:j+1], MixtureKickVarTail[:j+1],
                    dz[:j+1], rng, n_samples=N_GEOMETRY))
    for j in range(n_model):
        remaining_range = R0 - z[j]
        V_j, dz_j = V[:j+1], dz[:j+1]
        weights = mcs.weighted_eigenvalues(V_j, dz_j)
        probability = float(mcs.gchi2_exact_cdf(remaining_range, weights))
        reach_analytical[j] = probability
        reach_rng[j], reach_rng_error[j] = mcs.gchi2_cdf_rng(remaining_range, weights, rng, n_samples=N_RNG)
        reach_geometry[j], reach_geometry_error[j] = mcs.reach_probability_rng_geometry(remaining_range, V_j, dz_j, rng, n_samples=N_GEOMETRY)

    z_plot = np.r_[z, z[-1] + dz[-1]]
    dz_plot = np.r_[dz, dz[-1]]

    reach_g4_plot = np.r_[reach_analytical_g4, 0.0]
    reach_analytical_plot = np.r_[reach_analytical, 0.0]
    reach_rng_plot = np.r_[reach_rng, 0.0]
    reach_geometry_plot = np.r_[reach_geometry, 0.0]

    reach_rng_error_plot = np.r_[reach_rng_error, 0.0]
    reach_geometry_error_plot = np.r_[reach_geometry_error, 0.0]

    mixture_reach_analytical_plot = np.r_[mixture_reach_analytical, 0.0]
    mixture_reach_geometry_plot = np.r_[mixture_reach_geometry, 0.0]

    mixture_reach_geometry_error_plot = np.r_[mixture_reach_geometry_error, 0.0]
    
    print(f"Last model evaluation: {z[n_model-1]:.8f} cm; "
          f"gap to R0: {R0-z[n_model-1]:.8f} cm; "
          f"reach there: {reach_analytical[n_model-1]:.6g}")
    print(f"Geant4 reach at last recorded boundary: {reach_analytical_g4[-1]:.6g}")

    fig, (ax, ax_diff) = plt.subplots(
        2, 1, figsize=(12, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax.set_title("Reach probability vs. depth")
    ax.plot(z_plot, reach_g4_plot, ".-", color="black", label="Geant4")
    ax.plot(z_plot, reach_analytical_plot, ".-", color="blue", label="Single Gaussian CDF")
    
    
    if  ANGULAR_MODEL == "mixture":
        # ax.plot(z_plot, reach_geometry_plot, "-.", color="orange", label="Single Gaussian RNG")
        line, = ax.plot(z_plot, mixture_reach_analytical_plot, "-^", color="red",
                        label=f"Gaussian Mixture CDF ({N_CONFIGURATIONS} sampled configurations)")
        # ax.fill_between(z[:n_model],
        #                 np.maximum(0, mixture_reach_analytical[:n_model]
        #                            - 2*mixture_reach_analytical_error[:n_model]),
        #                 np.minimum(1, mixture_reach_analytical[:n_model]
        #                            + 2*mixture_reach_analytical_error[:n_model]),
        #                 color=line.get_color(), alpha=.2,
        #                 label="Configuration sampling ±2 SE (pointwise)")
        
        # ax.plot(z_plot, mixture_reach_geometry_plot, "-^", color="lime", label="Gaussian Mixture RNG")
    else:
        # ax.plot(z_plot, reach_rng_plot, ":.", color="blue", label="Reach Prob. - Chi-squared RNG")
        ax.plot(z_plot, reach_geometry_plot, "--o", color="orange", label="Reach Prob. - RNG")
        # ax.fill_between(z_plot, np.maximum(0, reach_rng_plot-2*reach_rng_error_plot),
                            # np.minimum(1, reach_rng_plot+2*reach_rng_error_plot), alpha=.25,
                            # label="RNG ±2 standard errors")
    ax.axvline(R0, color="gray", linestyle="--", label="R0")
    ax.set_ylabel("Reach probability")
    ax.legend()
    
    # ax.set_xscale("log")
    ax.grid()
    
    residuals = reach_analytical_g4 - reach_analytical
    mixture_residuals = reach_analytical_g4 - mixture_reach_analytical
            
    ax_diff.plot(z, residuals, "-o", color="blue", label="Geant4 - Single Gaussian CDF")
    # ax_diff.plot(z_plot, reach_geometry_plot-reach_analytical_plot, "-o", color="orange", label="RNG - Analytical")
    # ax_diff.fill_between(z_plot, -2*reach_rng_error_plot, 2*reach_rng_error_plot, alpha=.25, label="RNG - Analytical error band")
    if ANGULAR_MODEL == "mixture":
        ax_diff.plot(z, mixture_residuals, "-^", color="red", label="Geant4 - Gaussian Mixture sampled CDF")
        # ax_diff.plot(z_plot, mixture_reach_geometry_plot-mixture_reach_analytical_plot, "-^", color="lime", label="Gaussian Mixture RNG - Sampled CDF")
    ax_diff.axvline(R0, color="gray", linestyle="--", label="R0")
    ax_diff.set(xlabel="Depth / cm", ylabel="Residuals")
    ax_diff.grid()
    ax_diff.legend(loc='center left')
    fig.tight_layout()
    plt.xlim(left=X_LIM, right=z_plot[-1]+0.03)
    plt.savefig("MCSReachProbabilityPWO_mixture.svg", bbox_inches="tight")
    plt.show()

    # Densities live on the intervals between measured boundaries. Use these
    # exact same z edges; no centre grid, extra G4 edge, or normalization.
    # The first interval and the unobserved tail are deliberately not inferred.
    stop_bin_probabilityg4 = -np.diff(reach_g4_plot)
    stop_bin_densityg4 = stop_bin_probabilityg4 / dz_plot[1:] / 100
    
    stop_bin_probability = -np.diff(reach_analytical_plot)
    stop_bin_density = stop_bin_probability / dz_plot[1:] / 100
    
    stop_bin_probability_mixture = -np.diff(mixture_reach_analytical_plot)
    stop_bin_density_mixture = stop_bin_probability_mixture / dz_plot[1:] / 100
    
    if np.any(stop_bin_probability < -1e-8):
        warnings.warn("Analytical reach rises with depth; inspect CDF numerical stability.")
    plt.figure(figsize=(12, 9))
    if  ANGULAR_MODEL == "mixture":
        plt.stairs(np.r_[0, stop_bin_density_mixture, 0], np.r_[0, z_plot, 0], baseline=None, label="Mixture stopping probability")
    else:
        plt.stairs(np.r_[0, stop_bin_density, 0], np.r_[0, z_plot, 0], baseline=None, label="Analytical stopping probability")
    #plt.stairs(stop_analytical_plot, z_plot, baseline=None, label="Analytical PDF")
    plt.stairs(np.r_[0, stop_bin_densityg4, 0], np.r_[0, z_plot, 0], baseline=None, label="Geant4 stopping probability")
    plt.xlabel("Depth / cm")
    plt.ylabel(r"Stop Probability")# / cm$^{-1}$")
    plt.legend()
    plt.grid()
    plt.xlim(left=X_LIM, right=z_plot[-1]+0.03)
    plt.tight_layout()
    plt.savefig("MCSStopProbabilityPWO_mixture.svg", bbox_inches="tight")
    plt.show()

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

    plt.figure(figsize=(12, 9))
    plt.plot(z, np.degrees(np.sqrt(CumVarHighland)), label="Integrated Highland, global log")
    plt.scatter(z, np.degrees(CumSigmaCore), s=10, label="Geant4 cumulative core width")
    if ANGULAR_MODEL == "mixture":
        plt.plot(z, np.degrees(CumMixtureSigmaCore), label="Mixture narrow sigma")
        plt.plot(z, np.degrees(CumMixtureSigmaTail), label="Mixture broad sigma")
        plt.plot(z, np.degrees(np.sqrt(CumMixtureVariance)), "--", label="Mixture RMS")
    plt.xlabel("Depth / cm")
    plt.ylabel("Projected angle / degree")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("multiple_coulomb_scattering.svg", bbox_inches="tight")
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
    plt.savefig("lateral_scattering.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
