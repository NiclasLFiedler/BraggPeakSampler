import uproot
import numpy as np
#import matplotlib
#matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import gamma as gamma_func
from scipy.stats import chi2, gamma
from scipy.integrate import quad

from dataclasses import dataclass

import sys
sys.path.append("../../range_energy/data_analysis")
import analysisFunctions
import mcs_helper as mcs

@dataclass
class TargetParameters:
    Thickness: float
    Pmod_theo: float
    range: float
    sigma: float
    sigma_t: float
    t: float
    resRange: float
    sigma_res: float
    sigma_t_res: float
    t_res: float
    Pmod_res: float
    Pmod_sim: float
    energy: float
    sigma_E: float
    sigma_T_E: float
    t_E: float
    Pmod_E: float


def calculate_coefficients(weights):
    weights = np.asarray(weights, dtype=float)

    coefficients = np.ones(len(weights))

    for i, wi in enumerate(weights):
        for j, wj in enumerate(weights):
            if i != j:
                coefficients[i] *= wi / (wi - wj)

    return coefficients

def pdf_chi2_scaled(x, w):
    """PDF of scaled chi²_2: w * chi²_2"""
    return chi2.pdf(x / w, df=2) / w
 
def gchi2_exact_pdf(x, weights):
    weights = np.asarray(weights, dtype=float)
    weights = weights[weights > 0]

    coefficients = calculate_coefficients(weights)

    x = np.asarray(x, dtype=float)
    pdf = np.zeros_like(x)

    for wi, Ai in zip(weights, coefficients):
        pdf += Ai / (2 * wi) * np.exp(-x / (2 * wi))

    return pdf

def gchi2_exact_cdf(x, weights):
    weights = np.asarray(weights, dtype=float)
    weights = weights[weights > 0]

    coefficients = calculate_coefficients(weights)

    x = np.asarray(x, dtype=float)

    cdf = np.ones_like(x)

    for wi, Ai in zip(weights, coefficients):
        cdf -= Ai * np.exp(-x / (2 * wi))

    return cdf

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def gaussian_sigma_vs_depth(depth, angles, depths, bins=100):
    sigma_fit = np.full(len(depths), np.nan)
    variance_fit = np.full(len(depths), np.nan)

    std_data = np.full(len(depths), np.nan)
    variance_data = np.full(len(depths), np.nan)

    for i, d in enumerate(depths):

        selected = angles[np.abs(depth - d) < 0.001]
        selected = selected[np.isfinite(selected)]

        if len(selected) < 10:
            continue

        counts, bin_edges = np.histogram(selected, bins=bins, density=True)
        bin_centers = (0.5 * (bin_edges[:-1] + bin_edges[1:]))

        A0 = np.max(counts)
        mu0 = np.mean(selected)
        sigma0 = np.std(selected)

        try:
            popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=[A0, mu0, sigma0], maxfev=10000)
            A, mu, sigma = popt
            sigma = abs(sigma)
            sigma_fit[i] = sigma
            variance_fit[i] = sigma**2

        except RuntimeError:
            continue

        std_data[i] = np.std(selected)
        variance_data[i] = np.var(selected)

    return (sigma_fit, variance_fit, std_data, variance_data)

def plotSingleThickness(target_depth, depth, deltaX, label):
    selected_deltaX = deltaX[np.abs(depth - target_depth) < 0.001]
    selected_deltaX = selected_deltaX[np.isfinite(selected_deltaX)]

    plt.hist(selected_deltaX, bins=1000, density=True, alpha=0.7, label=label)

def gaussian_core_sigma(data):
    data = data[np.isfinite(data)]

    if len(data) < 10:
        return np.nan

    q_low, q_high = np.percentile( data, [15.8655, 84.1345])

    return 0.5 * (q_high - q_low)

with uproot.open("h2oproj.root") as f:
    tree = f["braggsampler"]

    event = tree["event"].array(library="np").astype(np.float32)
    layerID = tree["layerID"].array(library="np").astype(np.float32)
    depths = tree["depth"].array(library="np").astype(np.float32)
    CumScatteringAngle =   tree["CumScatteringAngle"].array(library="np").astype(np.float32)
    ScatteringAngle =   tree["SingleScatteringAngle"].array(library="np").astype(np.float32)
    deltaX =        tree["deltaX"].array(library="np").astype(np.float32)


event_ids, event_index = np.unique(event, return_inverse=True)

layers_int = layerID.astype(np.int64)
if not np.all(layerID == layers_int):
    raise ValueError("layerID contains noninteger values.")

max_layer = np.full(len(event_ids), -1, dtype=np.int64)
np.maximum.at(max_layer, event_index, layers_int)

last_row = np.r_[event[1:] != event[:-1], True]

print("Unique recorded events:", len(event_ids))
print("Contiguous event blocks:", np.count_nonzero(last_row))
print("Blocks ending below their event's maximum layer:",  np.count_nonzero(layers_int[last_row] < max_layer[event_index[last_row]])
)
del event_index
del event

z = np.unique(depths)
dz = np.diff(np.r_[0.0, z])
print("Layer thickness:", dz)

boundary_layers = np.arange(max_layer.max() + 2)
g4_reach_depths = (boundary_layers + 1) * dz
max_layer_sorted = np.sort(max_layer)

n_reaching = (len(max_layer_sorted) - np.searchsorted( max_layer_sorted, boundary_layers, side="left"))
reach_analytical_g4 = n_reaching / len(event_ids)

del max_layer_sorted
del event_ids

stop_bin_probabilityg4 = (reach_analytical_g4[:-1] - reach_analytical_g4[1:])
stop_bin_densityg4 = (stop_bin_probabilityg4 / np.diff(g4_reach_depths))
stop_bin_centresg4 = (g4_reach_depths[:-1] + g4_reach_depths[1:]) / 2

grad_stopping_probability_g4 = -np.gradient(reach_analytical_g4, g4_reach_depths)
grad_stopping_probability_g4 = np.maximum(grad_stopping_probability_g4, 0)


plt.figure(figsize=(10, 6))
plt.plot(g4_reach_depths, reach_analytical_g4, label="G4 Reach probability")
plt.plot(g4_reach_depths, grad_stopping_probability_g4, color="red", linewidth=2, label="G4 Stopping probability")
plt.plot(stop_bin_centresg4, stop_bin_probabilityg4, color="orange", linewidth=2, label="G4 Stopping bin differences")
plt.plot(stop_bin_centresg4, stop_bin_densityg4, color="green", linewidth=2, label="G4 stopping density — gradient")

plt.xlabel("Stopping depth / cm")
plt.ylabel(r"$P_{\mathrm{stop}}(x)$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

CumSigmaCore = np.full(len(z), np.nan)
SingleSigmaCore = np.full(len(z), np.nan)
CumDeltaXSigmaCore = np.full(len(z), np.nan)

print("Calculating Gaussian core sigma for each depth")

for i, d in enumerate(z):
    if i%10==0: print(f"Depth {d} cm of {z[-1]} cm")
    mask = np.abs(depths - d) < 0.001

    CumSigmaCore[i] = gaussian_core_sigma(CumScatteringAngle[mask])
    SingleSigmaCore[i] = gaussian_core_sigma(ScatteringAngle[mask])
    CumDeltaXSigmaCore[i] = gaussian_core_sigma(deltaX[mask])

del CumScatteringAngle
del ScatteringAngle
del deltaX
del depths

CumVarCore = CumSigmaCore**2
SingleAngleVarianceFromCum = np.empty_like(CumVarCore)
SingleAngleVarianceFromCum[0] = CumVarCore[0]
SingleAngleVarianceFromCum[1:] = np.maximum(CumVarCore[1:] - CumVarCore[:-1], 0)
SingleAngleRMSFromCum = np.sqrt(SingleAngleVarianceFromCum)

# ============================================================================
print("Calculating Depth-dependent generalized chi-square distribution")
# ============================================================================
usePbWO4 = False

data = analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/h2o_alt_range_energy.npz") if not usePbWO4 else analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/pbwo4_alt_range_energy.npz")

E0 = 220 
p_exp = data.p
alpha = data.alpha[0]

R0 = analysisFunctions.range_energy(data, E0)
R0 = 30.719

R0 = 30.75 if not usePbWO4 else 6.73398
X0 = 36.08 if not usePbWO4 else 0.89

z_mid = z - 0.5 * z

useMask = True
useHighland = False

if useMask:
    N = len(CumSigmaCore)
    z = z[:N]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[:N]
    SingleSigmaCore = SingleSigmaCore [:N]
    CumDeltaXSigmaCore = CumDeltaXSigmaCore [:N]
    
    # mask = z < R0 + 1
    
    # z = z[mask]
    # CumSigmaCore = CumSigmaCore[mask]
    # SingleAngleRMSFromCum = SingleAngleRMSFromCum[mask]
    # SingleSigmaCore = SingleSigmaCore[mask]
    # CumDeltaXSigmaCore = CumDeltaXSigmaCore[mask]


def energy_at_depth(depth_mid, E0, R0, p_exp):
    return E0 * (1.0 - depth_mid / R0)**(1.0 / p_exp)

if useHighland:
    CumVarHighland = mcs.highland_variance(z=depths, dz=dz, R0=R0, E0=E0, p_exp=p_exp, X0=X0, energy_at_depth= lambda mid: energy_at_depth(mid, E0, R0, p_exp))
    CumRMSHighland = np.sqrt(CumVarHighland)

z, dz, V = mcs.simulation_grid(z, CumVarCore, R0)

z, dz = mcs.two_grid(R0, R0-1.0, 0.10, 0.005)

rng = np.random.default_rng(12345)
reach_analytical = np.zeros(len(depthsFiltered))
reach_rng = np.zeros(len(depthsFiltered))
reach_rng_error = np.zeros(len(depthsFiltered))
reach_geometry = np.zeros(len(depthsFiltered))
reach_geometry_error = np.zeros(len(depthsFiltered))


if useHighland:
    print(
        f"Last calculated depth: {depthsFiltered[-1]:.6f} cm\n"
        f"R0: {R0:.6f} cm\n"
        f"Unresolved final interval: "
        f"{R0 - depthsFiltered[-1]:.6f} cm\n"
        f"Last reach probability: {reach_analytical[-1]:.6f}"
    )

    depthsFiltered = np.r_[0.0, depthsFiltered, R0]
    reach_analytical = np.r_[1.0, reach_analytical, 0.0]
    reach_rng = np.r_[1.0, reach_rng, 0.0]
    reach_rng_error = np.r_[0.0, reach_rng_error, 0.0]


if useHighland:
    depthsFiltered = depths.copy()
else:
    depthsFiltered = G4depths[valid_mask]

for j, d in enumerate(z):
        V_j, dz_j = V[:j+1], dz[:j+1]
        weights = weighted_eigenvalues(V_j, dz_j)
        # Your existing analytical CDF, not supplied by this module:
        reach_analytical[j] = (gchi2_exact_cdf(R0-d, weights)
                               if weights.size else 1.0)
        reach_rng[j], reach_rng_error[j] = gchi2_cdf_rng(
            R0-d, weights, rng)
        reach_geometry[j], reach_geometry_error[j], _ = (
            reach_probability_rng_geometry(R0-d, V_j, dz_j, rng))


for j, d in enumerate(depthsFiltered):
    deltaR_max = R0 - d

    if deltaR_max <= 0:
        continue

    V = CumVarHighland[:j+1] if useHighland else CumVarRad_filtered[:j+1]

    C_j = np.minimum.outer(V, V)
    eigenvalues = np.linalg.eigvalsh(C_j)
    eigenvalues = eigenvalues[eigenvalues > 0]

    integration_step = deltaz if useHighland else layerThickness

    weights = eigenvalues * integration_step / 2

    reach_analytical[j] = gchi2_exact_cdf(deltaR_max, weights)

    reach_rng[j], reach_rng_error[j] = mcs.gchi2_cdf_rng( deltaR_max, weights, rng, n_samples=100_000)

    reach_geometry[j], reach_geometry_error[j], _ = (mcs.reach_probability_rng_geometry( remaining_range=R0 - d, cumulative_variance=V, dz=deltaz, rng=rng, n_samples=50_000))

stopping_probability = -np.gradient(reach_analytical, depthsFiltered)
stopping_probability = np.maximum(stopping_probability, 0)

ig, (ax, ax_diff) = plt.subplots( 2, 1, figsize=(10, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
ax.plot( depthsFiltered, reach_analytical, label="Analytical CDF", linewidth=2)
ax.plot( depthsFiltered, reach_rng, "--", label="RNG CDF", linewidth=2)
ax.plot( depthsFiltered, reach_geometry, "--", label="RNG CDF without approximations", linewidth=2)
ax.fill_between( depthsFiltered, np.maximum(0, reach_rng - 2 * reach_rng_error), np.minimum(1, reach_rng + 2 * reach_rng_error), alpha=0.25, label="RNG ±2 standard errors")
ax.set_ylabel("Reach probability")
ax.legend()
ax.grid()

difference = reach_rng - reach_analytical
ax_diff.plot(depthsFiltered, difference, color="black")
ax_diff.fill_between( depthsFiltered, -2 * reach_rng_error, 2 * reach_rng_error, alpha=0.25)

ax_diff.axhline(0, color="gray", linestyle="--")
ax_diff.set_xlabel("Depth / cm")
ax_diff.set_ylabel("RNG - analytical")
ax_diff.grid()

plt.tight_layout()
plt.show()

##########################

plt.figure(figsize=(10, 6))
plt.plot(depthsFiltered, reach_analytical, linewidth=2, label="Analytical reach probability")
plt.plot(depthsFiltered, stopping_probability, linewidth=2, label="Analytical stopping probability")
plt.plot(depthsFiltered, reach_rng, "--", linewidth=2, label="RNG reach probability")
plt.plot(depthsFiltered, reach_geometry, "--", label="RNG CDF without approximations", linewidth=2)
plt.plot(g4_reach_depths, reach_analytical_g4, label="G4 Reach probability")
plt.plot(g4_reach_depths, grad_stopping_probability_g4, color="red", linewidth=2, label="G4 Stopping probability")
plt.axvline(R0, linestyle="--", color="black")
plt.xlabel("Stopping depth / cm")
plt.ylabel(r"$P_{\mathrm{stop}}(x)$")

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

#########################

plt.figure(figsize=(12, 9))

plt.plot(G4depths, SingleAngleRMSFromCum, "o-", label="Gaussian fit single scattering rms")
plt.plot(G4depths, SingleSigmaCore, "s--", label=r"Geant4")

plt.xlabel("Depth / cm")
plt.ylabel("Single scattering angle RMS / degree")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

################################################################################# gaussian test 

plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.plot(depths, np.degrees(np.sqrt(CumVarHighland)), color="navy", linewidth=2, label="Integral Highland global log")
plt.scatter(G4depths, CumSigmaCore, s=10, color="green", label="Geant4 Theta")
plt.scatter(G4depths, SingleAngleRMSFromCum, s=10, label="SingleAngleRMSFromCum")
plt.scatter(G4depths, SingleSigmaCore, s=10, label="SingleSigmaCore")

plt.xlabel("Depth / cm")
plt.ylabel("RMS projected angle / degree")
plt.title("Multiple Coulomb Scattering")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("multiple_coulomb_scattering.svg", format="svg", bbox_inches="tight")
plt.show()

layerThickness  = G4depths[1] - G4depths[0]
lateralVariance = np.zeros(len(G4depths))

for j, d in enumerate(G4depths):
    lever_arm = layerThickness+d - G4depths[:j+1]
    lateralVariance[j] = np.maximum(np.sum((lever_arm * np.tan(np.radians(SingleAngleRMSFromCum[:j+1])))**2),0)

plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.scatter(G4depths, np.sqrt(lateralVariance), marker="o", s=10, color="orange", label="Geant4")
plt.plot(G4depths, CumDeltaXSigmaCore, color="navy", linewidth=2, label="Geant4 Lateral Scattering RMS")
plt.xlabel("Depth / cm")
plt.ylabel("Lateral Scattering / cm")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("lateral_scattering.svg", format="svg", bbox_inches="tight")
plt.show()