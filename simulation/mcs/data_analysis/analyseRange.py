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



def calculate_eigenvalue_sum(eigenvalues):
    eigenvalues = np.asarray(eigenvalues, dtype=float)
    n = len(eigenvalues)
    result = np.zeros(n)
    
    for i in range(n):
        total = 0.0
        for j in range(n):
            if i != j:
                total += eigenvalues[i] / (eigenvalues[i] - eigenvalues[j])
        result[i] = total
    
    return result
 
 
def calculate_eigenvalue_sum_vectorized(eigenvalues):
    eigenvalues = np.asarray(eigenvalues, dtype=float)
    n = len(eigenvalues)

    numerator = eigenvalues[:, None]
    denominator = eigenvalues[:, None] - eigenvalues[None, :]
    
    with np.errstate(divide='ignore', invalid='ignore'):
        fraction_matrix = numerator / denominator
    
    np.fill_diagonal(fraction_matrix, 0)
    
    result = np.sum(fraction_matrix, axis=1)
    
    return result


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
    beta = w_sum_sq / w_sum
    
    return gamma.pdf(x, a=alpha, scale=beta)
 
# ============================================================================
# Method 3: Exact (via successive convolution)
# ============================================================================
def gchi2_exact_convolution(x, weights):
    """PDF using exact convolution of scaled chi²_2 distributions"""
    if len(weights) == 1:
        return pdf_chi2_scaled(x, weights[0])
    
    # Start with first chi²_2
    print(f"Convolving PDFs with number of weights: {len(weights)}")
    pdf_result = lambda t, w=weights[0]: pdf_chi2_scaled(t, w)
    
    # Convolve with each subsequent chi²_2
    for i, w in enumerate(weights[1:]):
        pdf_prev = pdf_result
        print(f"Convolving with weight: i={i+1}")
        pdf_result = lambda t, w=w, prev=pdf_prev: conv_two_pdfs(prev, pdf_chi2_scaled, t, 1, w)
    
    return pdf_result(x)


def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def gaussian_sigma_vs_depth(depth, angles, depths, tolerance=0.1, bins=100):
    sigma_fit = np.full(len(depths), np.nan)
    variance_fit = np.full(len(depths), np.nan)

    std_data = np.full(len(depths), np.nan)
    variance_data = np.full(len(depths), np.nan)

    for i, d in enumerate(depths):

        selected = angles[np.abs(depth - d) < tolerance]
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
    tolerance = 0.1

    selected_deltaX = deltaX[np.abs(depth - target_depth) < tolerance]
    selected_deltaX = selected_deltaX[np.isfinite(selected_deltaX)]

    plt.hist(selected_deltaX, bins=1000, density=True, alpha=0.7, label=label)

with uproot.open("h2oproj.root") as f:
    tree = f["braggsampler"]

    event = tree["event"].array(library="np")
    depth = tree["depth"].array(library="np")
    CumScatteringAngle =   np.degrees(tree["CumScatteringAngle"].array(library="np"))
    ScatteringAngle =   np.degrees(tree["SingleScatteringAngle"].array(library="np"))
    deltaX =        tree["deltaX"].array(library="np")

G4depths = np.unique(depth)

CumSigmaFit, CumVarFit, CumStd, _ =  gaussian_sigma_vs_depth( depth, CumScatteringAngle, G4depths, tolerance=0.1, bins=2000)

SingleSigmaFit, SingleVarFit, SingleStd, _ =  gaussian_sigma_vs_depth(depth, ScatteringAngle, G4depths, tolerance=0.1, bins=2000)

CumDeltaXSigma, CumDeltaXVar, CumDeltaXStd, _ =  gaussian_sigma_vs_depth(depth, deltaX, G4depths, tolerance=0.1, bins=2000)

SingleAngleVarianceFromCum = np.empty_like(CumVarFit)
SingleAngleVarianceFromCum[0] = CumVarFit[0]
SingleAngleVarianceFromCum[1:] = np.maximum(CumVarFit[1:] - CumVarFit[:-1], 0)
SingleAngleRMSFromCum = np.sqrt(SingleAngleVarianceFromCum)

plt.figure(figsize=(12, 9))

plt.plot(G4depths, SingleAngleRMSFromCum, "o-", label="Gaussian fit single scattering rms")
plt.plot(G4depths, SingleSigmaFit, "s--", label=r"$\sqrt{\mathrm{np.var}}$")

plt.xlabel("Depth / mm")
plt.ylabel("Single scattering angle RMS / degree")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

target_depth = 20.0
tolerance = 0.1

angles = ScatteringAngle[np.abs(depth - target_depth) < tolerance]
counts, bin_edges = np.histogram(angles, bins=500, density=True)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
A0 = np.max(counts)
mu0 = np.mean(angles)
sigma0 = np.std(angles)
popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=[A0, mu0, sigma0], maxfev=100000)
A, mu, sigma = popt
x_fit = np.linspace( bin_edges[0], bin_edges[-1], 500)
y_fit = gaussian(x_fit, A, mu, sigma)

data = analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/h2o_alt_range_energy.npz")
alpha = data.alpha[0]
p_exp = data.p
E0 = 220

R0 = analysisFunctions.range_energy(data, E0)

useMask = True
if useMask:
    N = len(CumSigmaFit)
    G4depths = G4depths[:N]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[:N]
    SingleSigmaFit = SingleSigmaFit [:N]
    CumDeltaXSigma = CumDeltaXSigma [:N]
    
    mask = G4depths < R0 + 1

    G4depths = G4depths[mask]
    CumSigmaFit = CumSigmaFit[mask]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[mask]
    SingleSigmaFit = SingleSigmaFit[mask]
    CumDeltaXSigma = CumDeltaXSigma[mask]

m_p = 938.272
X0 = 36.08
dx = 1
x_max = 0.999 * R0
depths = np.arange(dx, x_max, dx)

E_k = E0 * (1 - depths/R0)**(1/p_exp)
betaPc = E_k * (E_k + 2*m_p)/(E_k + m_p)

integrand = (13.6/betaPc)**2 * dx/X0
log_factor = 1 + 0.038*np.log(depths/X0)
theta_integrated_rad = log_factor * np.sqrt(np.cumsum(integrand))
theta_integrated_deg = np.degrees(theta_integrated_rad)

varianceHigh = theta_integrated_rad**2
SingleAngleVarianceHigh = np.empty_like(theta_integrated_deg)
SingleAngleVarianceHigh[0] = varianceHigh[0]
SingleAngleVarianceHigh[1:] = (varianceHigh[1:] - varianceHigh[:-1])
SingleAngleHigh = np.degrees(np.sqrt(SingleAngleVarianceHigh))

# ################
# integrand14 = (14.1/betaPc)**2 * dx/X0
# log_factor14 = 1 + 1/9*np.log(depths/X0)
# theta_integrated_deg14 = np.degrees(log_factor14 * np.sqrt(np.cumsum(integrand14)))

# pv0 = betaPc[0]       # approximately initial pv
# pv  = betaPc
# f_dM = (
    # 0.5244
    # + 0.1975*np.log10(1 - (pv/pv0)**2)
    # + 0.2320*np.log10(pv)
    # - 0.0098*np.log10(pv)
    #   *np.log10(1 - (pv/pv0)**2)
# )

# T = f_dM * (E_k/pv)**2 / X0

# variance = np.cumsum(T * dx)
# theta = np.degrees(np.sqrt(variance))

# theta_step = (13.6 / betaPc) * np.sqrt(dx / X0) * (1 + 0.038 * np.log(dx / X0))
# theta_naive_deg = np.degrees(np.sqrt(np.cumsum(theta_step**2)))

# simpleTheta = np.degrees((13.6 / betaPc) * np.sqrt(depths / X0) * (1 + 0.038 * np.log(depths / X0)))
# #############

plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.plot(depths, theta_integrated_deg, color="navy", linewidth=2, label="Integral Highland global log")

# plt.plot(depths, theta_integrated_deg14, color="green", linewidth=2, label="Integral Highland global log 14.1")
#plt.plot(depths, theta_naive_deg, color="orange", linewidth=2, label="Integral Highland local log")
# plt.plot(depths, theta, color="red", linewidth=2, label="Simple Highland")
#plt.plot(depths, SingleAngleHigh, color="black", linewidth=2.5, label="Single Highland Angle")

plt.scatter(G4depths, CumSigmaFit, s=10, color="green", label="Geant4 Theta")
plt.scatter(G4depths, SingleAngleRMSFromCum, s=10, label="SingleAngleRMSFromCum")
plt.scatter(G4depths, SingleSigmaFit, s=10, label="SingleSigmaFit")

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
plt.plot(G4depths, CumDeltaXSigma, color="navy", linewidth=2, label="Geant4 Lateral Scattering RMS")
plt.xlabel("Depth / cm")
plt.ylabel("Lateral Scattering / cm")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("lateral_scattering.svg", format="svg", bbox_inches="tight")
plt.show()

rangeDecrease = layerThickness/np.cos(np.radians(CumScatteringAngle))-layerThickness
rangeDecreaseSigma = np.sqrt(np.var(rangeDecrease))
sigmaRangeHighland = layerThickness/np.sqrt(2)*np.radians(theta_integrated_deg)**2
sigmaRange = layerThickness/np.sqrt(2)*np.radians(CumSigmaFit)**2

print(f"length range: {len(sigmaRange)}, length highland {len(sigmaRangeHighland)}")

plt.figure(figsize=(10, 7))
#plotSingleThickness(20, depth, rangeDecrease, r"$z/\cos(\theta)-z$")
plt.plot(G4depths, sigmaRange, color="navy", linestyle="--", marker="o", linewidth=2, label="Geant4 Lateral Scattering RMS")
plt.plot(depths, sigmaRangeHighland, color="red", linewidth=2, label="Highland Lateral Scattering RMS")
# plt.plot(G4depths, rangeDecreaseSigma, color="navy", linewidth=2, label="rangeDecreaseSigma")
plt.xlabel("Depth / cm")
plt.ylabel("RMS range decrease / cm")
# plt.yscale("log")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

valid_mask = ~np.isnan(CumVarFit)
CumVarFit_filtered = CumVarFit[valid_mask]
C = np.minimum.outer(CumVarFit_filtered, CumVarFit_filtered)

eigenvalues, eigenvectors = np.linalg.eigh(C)

depthsFiltered = depths[:len(CumVarFit_filtered)]
weights = eigenvalues*layerThickness/2
print(f"weights: {weights}")
a_is = [calculate_eigenvalue_sum(weights[:idx+1]) for idx in range(len(weights))]

x = depths
PDF_matrix_1 = np.zeros((len(depthsFiltered), len(weights)))

for idx, d in enumerate(depthsFiltered):
    a_i_truncated = a_is[idx]  # This has length idx+1
    for j in range(len(a_i_truncated)):
        PDF_matrix[idx, j] = a_i_truncated[j] * (1 - np.exp(-d / (2*weights[j])))


print(PDF_reach)
# a_i2 = calculate_eigenvalue_sum_vectorized(eigenvalues)
plt.figure(figsize=(10, 6))

plt.plot(depthsFiltered, PDF_reach, 'r--', linewidth=2)

plt.xlabel('x', fontsize=12)
plt.ylabel('PDF', fontsize=12)
plt.tight_layout()
plt.show()

pdf_sat = np.array([gchi2_satterthwaite(xi, weights) for xi in depthsFiltered])
print(f"Computed Satterthwaite PDF for {len(x)} points.")
pdf_ww = np.array([gchi2_welch_welford(xi, weights) for xi in depthsFiltered])
print(f"Computed Welch-Welford PDF for {len(x)} points.")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(depthsFiltered, pdf_sat, 'b-', label='Satterthwaite (χ²)', linewidth=2)
plt.plot(depthsFiltered, pdf_ww, 'g-', label='Welch-Welford (Gamma)', linewidth=2)
plt.plot(depthsFiltered, PDF_reach, 'r--', label='Exact (Convolution)', linewidth=2.5)

plt.xlabel('x', fontsize=12)
plt.ylabel('PDF', fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

