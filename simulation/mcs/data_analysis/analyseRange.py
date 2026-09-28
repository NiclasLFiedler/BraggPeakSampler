"""One Geant4 boundary grid for reach, scattering and lateral plots.

Requires your existing analysisFunctions.py, mcs_helper.py, ROOT input and
range calibration. All depth coordinates are cm; input angles are radians.
Scoring must identify actual downstream crossings, with one nominal depth
per layer and contiguous layer IDs starting at zero. No interpolation.
"""
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import uproot

sys.path.append("../../range_energy/data_analysis")
import analysisFunctions
import mcs_helper as mcs

# Configuration: preserve the settings in the supplied script.
usePbWO4 = False
useHighland = False
ROOT_FILE = "h2oproj.root"
E0 = 220.0
# Set to None to use range_energy(data, E0). These are your current overrides.
R0_OVERRIDE = 6.73398 if usePbWO4 else None #30.72
# Set to the generated-primary count to include events with no recorded hit.
# None preserves your original normalization to recorded events only.
N_PRIMARIES = None
N_RNG = 100_000
N_GEOMETRY = 100_000
# Confirm against your scorer: convert deltaX to cm if it was stored in mm.
DELTA_X_TO_CM = 1.0


def gaussian_core_sigma(values, latestSigma):
    values = values[np.isfinite(values)]
    if len(values) < 100:
        return latestSigma*1.01
    low, high = np.percentile(values, [15.8655, 84.1345])
    return 0.5 * (high - low)


def integer_ids(values, name):
    values = np.asarray(values)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError(f"{name} must contain finite nonnegative IDs.")
    ids = values.astype(np.int64)
    if not np.all(values == ids):
        raise ValueError(f"{name} contains noninteger values.")
    return ids

def main():
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
    
    print("Calculating angular core widths by layer")
    
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
        
        if i % 10 == 0:
            print(f"Layer {i}, depth {z[i]:.6f} cm")
            
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
                 if usePbWO4 else "../../range_energy/data_analysis/h2o_alt_range_energy.npz")
    data = analysisFunctions.load_EnergyRange(data_file)
    p_exp = float(data.p)
    fitted_R0 = float(analysisFunctions.range_energy(data, E0))
    R0 = fitted_R0 if R0_OVERRIDE is None else float(R0_OVERRIDE)
    X0 = 0.89 if usePbWO4 else 36.08
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
    CumVarHighland[:n_model] = mcs.highland_variance(
        z=z[:n_model], dz=dz[:n_model], R0=R0,
        E0=E0, p_exp=p_exp, X0=X0,
    )
    # Above retains your existing power-law energy profile. For your polynomial
    # law, pass the consistent inverse as energy_at_depth to the helper.
    V = CumVarHighland if useHighland else CumVarCore
    # Validate without replacing/truncating the master z/dz arrays.
    mcs.simulation_grid(z[:n_model], V[:n_model], R0,
                        layer_ids=boundary_layers[:n_model])

    rng = np.random.default_rng(12345)
    reach_analytical = np.zeros(len(z))
    stop_analytical = np.zeros(len(z))
    reach_rng = np.zeros(len(z))
    reach_rng_error = np.zeros(len(z))
    reach_geometry = np.zeros(len(z))
    reach_geometry_error = np.zeros(len(z))
    
    for j in range(n_model):
        remaining_range = R0 - z[j]
        V_j, dz_j = V[:j+1], dz[:j+1]
        weights = mcs.weighted_eigenvalues(V_j, dz_j)
        probability = float(mcs.gchi2_exact_cdf(remaining_range, weights))
        stoppProbability = float(mcs.gchi2_exact_pdf(remaining_range, weights))
        if not np.isfinite(probability) or not 0 <= probability <= 1:
            raise FloatingPointError(
                f"Analytical CDF={probability} at z={z[j]}; check residue cancellation.")
        reach_analytical[j] = probability
        stop_analytical[j] = stoppProbability
        reach_rng[j], reach_rng_error[j] = mcs.gchi2_cdf_rng(
            remaining_range, weights, rng, n_samples=N_RNG)
        # Indexing the first two results also accepts your two-return variant.
        geometry_result = mcs.reach_probability_rng_geometry(
            remaining_range, V_j, dz_j, rng, n_samples=N_GEOMETRY)
        reach_geometry[j], reach_geometry_error[j] = geometry_result

    z_plot = np.r_[z, z[-1] + dz[-1]]
    dz_plot = np.r_[dz, dz[-1]]

    reach_g4_plot = np.r_[reach_analytical_g4, 0.0]
    reach_analytical_plot = np.r_[reach_analytical, 0.0]
    stop_analytical_plot = np.r_[stop_analytical, 0.0]
    reach_rng_plot = np.r_[reach_rng, 0.0]
    reach_geometry_plot = np.r_[reach_geometry, 0.0]

    reach_rng_error_plot = np.r_[reach_rng_error, 0.0]
    reach_geometry_error_plot = np.r_[reach_geometry_error, 0.0]

    # Values at z >= R0 remain the fixed-range model's zero boundary values.
    print(f"Last model evaluation: {z[n_model-1]:.8f} cm; "
          f"gap to R0: {R0-z[n_model-1]:.8f} cm; "
          f"reach there: {reach_analytical[n_model-1]:.6g}")
    print(f"Geant4 reach at last recorded boundary: {reach_analytical_g4[-1]:.6g}")

    fig, (ax, ax_diff) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(z_plot, reach_analytical_plot, label="Analytical CDF")
    ax.plot(z_plot, reach_rng_plot, "--", label="Chi-squared RNG")
    ax.plot(z_plot, reach_geometry_plot, ":", label="Nonlinear-geometry RNG")
    ax.plot(z_plot, reach_g4_plot, ".-", label="Geant4 reach")
    ax.fill_between(z_plot, np.maximum(0, reach_rng_plot-2*reach_rng_error_plot),
                    np.minimum(1, reach_rng_plot+2*reach_rng_error_plot), alpha=.25,
                    label="RNG ±2 standard errors")
    ax.axvline(R0, color="gray", linestyle="--", label="R0")
    ax.set_ylabel("Reach probability")
    ax.legend()
    ax.grid()
    ax_diff.plot(z_plot, reach_rng_plot-reach_analytical_plot, color="black")
    ax_diff.fill_between(z_plot, -2*reach_rng_error_plot, 2*reach_rng_error_plot, alpha=.25)
    ax_diff.axhline(0, color="gray", linestyle="--")
    ax_diff.set(xlabel="Depth / cm", ylabel="RNG − analytical")
    ax_diff.grid()
    fig.tight_layout()
    plt.show()

    # Densities live on the intervals between measured boundaries. Use these
    # exact same z edges; no centre grid, extra G4 edge, or normalization.
    # The first interval and the unobserved tail are deliberately not inferred.
    stop_bin_probabilityg4 = -np.diff(reach_g4_plot)
    stop_bin_densityg4 = stop_bin_probabilityg4 / dz_plot[1:]
    stop_bin_probability = -np.diff(reach_analytical_plot)
    stop_bin_density = stop_bin_probability / dz_plot[1:]
    if np.any(stop_bin_probability < -1e-8):
        warnings.warn("Analytical reach rises with depth; inspect CDF numerical stability.")
    plt.figure(figsize=(10, 6))
    plt.stairs(stop_bin_density, z_plot, baseline=None, label="Analytical stopping density")
    #plt.stairs(stop_analytical_plot, z_plot, baseline=None, label="Analytical PDF")
    plt.stairs(stop_bin_densityg4, z_plot, baseline=None, label="G4 last-crossing interval density")
    plt.xlabel("Depth / cm")
    plt.ylabel(r"Density / cm$^{-1}$")
    plt.legend()
    plt.grid()
    plt.tight_layout()
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
    plt.scatter(z, np.degrees(SingleAngleRMSFromCum), s=10, label="From cumulative variance")
    plt.scatter(z, np.degrees(SingleSigmaCore), s=10, label="Geant4 single-angle core width")
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
