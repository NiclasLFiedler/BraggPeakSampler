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
 
def conv_two_pdfs(pdf1, pdf2, x, a, b):
    """Convolve two PDFs numerically"""
    def integrand(t):
        return pdf1(t, a) * pdf2(x - t, b)
    
    result, _ = quad(integrand, 0, max(x, 1e-10), limit=100)
    return result
 
# ============================================================================
# Method 1: Satterthwaite Approximation
# ============================================================================
def gchi2_satterthwaite(x, weights):
    """PDF using Satterthwaite approximation (scaled chi²)"""
    w_sum = np.sum(weights)
    w_sum_sq = np.sum(np.array(weights) ** 2)
    
    # Effective DOF and scale
    nu = 2 * w_sum**2 / w_sum_sq
    c = w_sum_sq / w_sum
    
    return chi2.pdf(x / c, df=nu) / c
 
# ============================================================================
# Method 2: Welch-Welford (Gamma Approximation)
# ============================================================================
def gchi2_welch_welford(x, weights):
    """PDF using Welch-Welford approximation (gamma distribution)"""
    w_sum = np.sum(weights)
    w_sum_sq = np.sum(np.array(weights) ** 2)
    
    # Gamma parameters
    alpha = w_sum**2 / w_sum_sq
    beta = 2*w_sum_sq / w_sum
    
    return gamma.pdf(x, a=alpha, scale=beta)
 
# ============================================================================
# Method 3: Exact (via successive convolution)
# ============================================================================
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

def gchi2_cdf_rng(x, weights, rng, n_samples=50_000, batch_size=1000):
    """Estimate P(sum_k weights[k] * chi2_2[k] <= x)."""
    weights = np.asarray(weights, dtype=float)

    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("Weights must be finite and nonnegative.")

    weights = weights[weights > 0]

    if x < 0:
        return 0.0, 0.0

    if weights.size == 0:
        return 1.0, 0.0

    n_reach = 0

    for start in range(0, n_samples, batch_size):
        size = min(batch_size, n_samples - start)

        # Independent chi² variables with 2 degrees of freedom
        samples = rng.chisquare(
            df=2,
            size=(size, len(weights))
        )

        delta_R = samples @ weights
        n_reach += np.count_nonzero(delta_R <= x)

    probability = n_reach / n_samples

    error = np.sqrt(probability * (1.0 - probability) / n_samples)

    return probability, error

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

    q_low, q_high = np.percentile(
        data,
        [15.8655, 84.1345]
    )

    return 0.5 * (q_high - q_low)

with uproot.open("h2oproj.root") as f:
    tree = f["braggsampler"]

    event = tree["event"].array(library="np")
    layerID = tree["layerID"].array(library="np")
    depth = tree["depth"].array(library="np")
    CumScatteringAngle =   np.degrees(tree["CumScatteringAngle"].array(library="np"))
    ScatteringAngle =   np.degrees(tree["SingleScatteringAngle"].array(library="np"))
    deltaX =        tree["deltaX"].array(library="np")


G4depths = np.unique(depth)
layerThickness  = G4depths[1] - G4depths[0]
print("Layer thickness:", layerThickness)

last = np.r_[event[1:] != event[:-1], True]
last_event = event[last]
last_layer = layerID[last]
last_depth = depth[last]
layers, counts = np.unique(last_layer, return_counts=True)
print("Number of total events:", len(event), " Number of last events:", len(last_event), "Number of layers:", len(last_layer), "Number of depths:", len(last_depth))

reach_analytical_g4 = np.array([np.mean(last_depth >= d) for d in G4depths])

stopping_probability_g4 = -np.gradient(reach_analytical_g4, G4depths)
stopping_probability_g4 = np.maximum(stopping_probability_g4, 0)
normalization = np.trapezoid(stopping_probability_g4, G4depths)

if normalization > 0:
    stopping_probability_g4 /= normalization

CumSigmaCore = np.full(len(G4depths), np.nan)
SingleSigmaCore = np.full(len(G4depths), np.nan)
CumDeltaXSigmaCore = np.full(len(G4depths), np.nan)

print("Calculating Gaussian core sigma for each depth")

# for i, d in enumerate(G4depths):
    # print(f"Depth {d} cm of {G4depths[-1]} cm")
    # selected = CumScatteringAngle[np.abs(depth - d) < 0.001]
    # CumSigmaCore[i] = gaussian_core_sigma(selected)

    # selected = ScatteringAngle[np.abs(depth - d) < 0.001]
    # SingleSigmaCore[i] = gaussian_core_sigma(selected)

    # selected = deltaX[np.abs(depth - d) < 0.001]
    # CumDeltaXSigmaCore[i] = gaussian_core_sigma(selected)

CumVarCore = CumSigmaCore**2

# ============================================================================
print("Calculating Depth-dependent generalized chi-square distribution")
# ============================================================================
data = analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/h2o_alt_range_energy.npz")
p_exp = data.p
alpha = data.alpha[0]

E0 = 220
R0 = analysisFunctions.range_energy(data, E0)

SingleAngleVarianceFromCum = np.empty_like(CumVarCore)
SingleAngleVarianceFromCum[0] = CumVarCore[0]
SingleAngleVarianceFromCum[1:] = np.maximum(CumVarCore[1:] - CumVarCore[:-1], 0)
SingleAngleRMSFromCum = np.sqrt(SingleAngleVarianceFromCum)

useMask = True
if useMask:
    N = len(CumSigmaCore)
    G4depths = G4depths[:N]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[:N]
    SingleSigmaCore = SingleSigmaCore [:N]
    CumDeltaXSigmaCore = CumDeltaXSigmaCore [:N]
    
    mask = G4depths < R0 + 1

    G4depths = G4depths[mask]
    CumSigmaCore = CumSigmaCore[mask]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[mask]
    SingleSigmaCore = SingleSigmaCore[mask]
    CumDeltaXSigmaCore = CumDeltaXSigmaCore[mask]

m_p = 938.272
X0 = 36.08
dx = layerThickness
x_max = 1 * R0
depths = np.arange(dx, x_max, dx)

E_k = E0 * (1 - depths/R0)**(1/p_exp)
betaPc = E_k * (E_k + 2*m_p)/(E_k + m_p)

integrand = (13.6/betaPc)**2 * dx/X0
log_factor = 1 + 0.038*np.log(depths/X0)
CumRMSHighland = log_factor * np.sqrt(np.cumsum(integrand))
CumVarHighland = CumRMSHighland**2
print(f"Projected Max Range: {R0:.3f} cm")

CumVarRad = np.radians(CumSigmaCore)**2
valid_mask = np.isfinite(CumVarRad)

depthsFiltered = G4depths[:len(CumVarRad)][valid_mask]
CumVarRad_filtered = CumVarRad[valid_mask]

layerThickness = G4depths[1] - G4depths[0]

useHighland = True
if useHighland:
    depthsFiltered = depths

rng = np.random.default_rng(12345)
reach_analytical = np.zeros(len(depthsFiltered))
reach_rng = np.zeros(len(depthsFiltered))
reach_rng_error = np.zeros(len(depthsFiltered))

storeEigenvalues = False

all_weights = []
all_eigenvalues = []


for j, d in enumerate(depthsFiltered):
    deltaR_max = R0 - d

    if deltaR_max <= 0:
        continue

    if useHighland:
        V = CumVarHighland[:j+1]
    else:
        V = CumVarRad_filtered[:j+1]

    C_j = np.minimum.outer(V, V)
    eigenvalues = np.linalg.eigvalsh(C_j)
    eigenvalues = eigenvalues[eigenvalues > 0]
    weights = eigenvalues * layerThickness / 2

    #mean_direct = layerThickness * np.sum(V)
    #mean_eigen = 2 * np.sum(weights)
    # print("mean_direct:", mean_direct, "mean_eigen:", mean_eigen)

    if storeEigenvalues:
        all_eigenvalues.append(eigenvalues)
        all_weights.append(weights)

    reach_analytical[j] = gchi2_exact_cdf(deltaR_max, weights)

    reach_rng[j], reach_rng_error[j] = gchi2_cdf_rng( deltaR_max, weights, rng, n_samples=100_000)


stopping_probability = -np.gradient(reach_analytical, depthsFiltered)
stopping_probability = np.maximum(stopping_probability, 0)
normalization = np.trapezoid(stopping_probability, depthsFiltered)
if normalization > 0:
    stopping_probability /= normalization

ig, (ax, ax_diff) = plt.subplots( 2, 1, figsize=(10, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})

ax.plot( depthsFiltered, reach_analytical, label="Analytical CDF", linewidth=2)
ax.plot( depthsFiltered, reach_rng, "--", label="RNG CDF", linewidth=2)
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
plt.plot(G4depths, reach_analytical_g4, label="G4 Reach probability")
plt.plot(G4depths, stopping_probability_g4, color="red", linewidth=2, label="G4 Stopping probability")

plt.xlabel("Stopping depth / cm")
plt.ylabel(r"$P_{\mathrm{stop}}(x)$")

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

#########################

plt.figure(figsize=(12, 9))

plt.plot(G4depths, SingleAngleRMSFromCum, "o-", label="Gaussian fit single scattering rms")
plt.plot(G4depths, SingleSigmaCore, "s--", label=r"$\sqrt{\mathrm{np.var}}$")

plt.xlabel("Depth / cm")
plt.ylabel("Single scattering angle RMS / degree")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

################################################################################# gaussian test 

# target_depth = 20.0
# angles = ScatteringAngle[np.abs(depth - target_depth) < 0.001]
# counts, bin_edges = np.histogram(angles, bins=500, density=True)
# bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
# A0 = np.max(counts)
# mu0 = np.mean(angles)
# sigma0 = np.std(angles)
# popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=[A0, mu0, sigma0], maxfev=100000)
# A, mu, sigma = popt
# x_fit = np.linspace( bin_edges[0], bin_edges[-1], 500)
# y_fit = gaussian(x_fit, A, mu, sigma)

alpha = data.alpha[0]
p_exp = data.p
E0 = 220

R0 = analysisFunctions.range_energy(data, E0)

useMask = True
if useMask:
    N = len(CumSigmaCore)
    G4depths = G4depths[:N]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[:N]
    SingleSigmaCore = SingleSigmaCore [:N]
    CumDeltaXSigmaCore = CumDeltaXSigmaCore [:N]
    
    mask = G4depths < R0 + 1

    G4depths = G4depths[mask]
    CumSigmaCore = CumSigmaCore[mask]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[mask]
    SingleSigmaCore = SingleSigmaCore[mask]
    CumDeltaXSigmaCore = CumDeltaXSigmaCore[mask]

m_p = 938.272
X0 = 36.08
dx = 1
x_max = 0.999 * R0
depths = np.arange(dx, x_max, dx)

E_k = E0 * (1 - depths/R0)**(1/p_exp)
betaPc = E_k * (E_k + 2*m_p)/(E_k + m_p)

integrand = (13.6/betaPc)**2 * dx/X0
log_factor = 1 + 0.038*np.log(depths/X0)
CumRMSHighland = log_factor * np.sqrt(np.cumsum(integrand))
theta_integrated_deg = np.degrees(CumRMSHighland)

varianceHigh = CumRMSHighland**2
SingleAngleVarianceHigh = np.empty_like(theta_integrated_deg)
SingleAngleVarianceHigh[0] = varianceHigh[0]
SingleAngleVarianceHigh[1:] = (varianceHigh[1:] - varianceHigh[:-1])
SingleAngleHigh = np.degrees(np.sqrt(SingleAngleVarianceHigh))

plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.plot(depths, theta_integrated_deg, color="navy", linewidth=2, label="Integral Highland global log")

# plt.plot(depths, theta_integrated_deg14, color="green", linewidth=2, label="Integral Highland global log 14.1")
#plt.plot(depths, theta_naive_deg, color="orange", linewidth=2, label="Integral Highland local log")
# plt.plot(depths, theta, color="red", linewidth=2, label="Simple Highland")
#plt.plot(depths, SingleAngleHigh, color="black", linewidth=2.5, label="Single Highland Angle")

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


pdf_sat = np.array([gchi2_satterthwaite(xi, weights) for xi in depthsFiltered])
pdf_ww = np.array([gchi2_welch_welford(xi, weights) for xi in depthsFiltered])

plt.figure(figsize=(10, 6))
plt.plot(depthsFiltered, pdf_sat, 'b-', label='Satterthwaite (χ²)', linewidth=2)
plt.plot(depthsFiltered, pdf_ww, 'g-', label='Welch-Welford (Gamma)', linewidth=2)

plt.xlabel('x', fontsize=12)
plt.ylabel('PDF', fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
