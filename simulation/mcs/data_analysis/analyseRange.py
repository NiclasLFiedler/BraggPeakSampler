import uproot
import numpy as np
#import matplotlib
#matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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
popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=[A0, mu0, sigma0])
A, mu, sigma = popt
x_fit = np.linspace( bin_edges[0], bin_edges[-1], 500)
y_fit = gaussian(x_fit, A, mu, sigma)

data = analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/h2o_alt_range_energy.npz")
alpha = data.alpha[0]
p_exp = data.p
E0 = 220

R0 = analysisFunctions.range_energy(data, E0)
print(f"range {R0}")
useMask = False
if useMask:
    N = len(CumSigmaFit)
    G4depths = G4depths[:N]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[:N]
    SingleSigmaFit = SingleSigmaFit [:N]

    mask = G4depths < R0 + 1

    G4depths = G4depths[mask]
    CumSigmaFit = CumSigmaFit[mask]
    SingleAngleRMSFromCum = SingleAngleRMSFromCum[mask]
    SingleSigmaFit = SingleSigmaFit[mask]

m_p = 938.272
X0 = 36.08
dx = 0.2
x_max = 0.999 * R0
depths = np.arange(dx, x_max, dx)

E_k = E0 * (1 - depths/R0)**(1/p_exp)
beta_p = E_k * (E_k + 2*m_p)/(E_k + m_p)

integrand = (13.6/beta_p)**2 * dx/X0
uncorrected_variance = np.cumsum(integrand)
log_factor = 1 + 0.038*np.log(depths/X0)
theta_integrated_rad = log_factor * np.sqrt(uncorrected_variance)
theta_integrated_deg = np.degrees(theta_integrated_rad)

varianceHigh = theta_integrated_rad**2
SingleAngleVarianceHigh = np.empty_like(theta_integrated_deg)
SingleAngleVarianceHigh[0] = varianceHigh[0]
SingleAngleVarianceHigh[1:] = (varianceHigh[1:] - varianceHigh[:-1])

SingleAngleHigh = np.degrees(np.sqrt(SingleAngleVarianceHigh))

log_step = 1 + 0.038 * np.log(dx / X0)
theta_step = (13.6 / beta_p) * np.sqrt(dx / X0) * log_step
theta_naive_deg = np.degrees(np.sqrt(np.cumsum(theta_step**2)))

print(f"Depth of material: {depths[-1]:.4f} cm")
print(f"Final cumulative scattering angle: {theta_integrated_deg[-1]:.4f} mrad")

plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.plot(depths, theta_integrated_deg, color="navy", linewidth=2, label="Integral Highland (Thick Target)")
plt.plot(depths, SingleAngleHigh, color="black", linewidth=2.5, label="Single Highland Angle")

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
#    print(f"Depth: {d:.2f} cm, Lever arm: {lever_arm}, lateralVariance[j]: {lateralVariance[j]:.4f} cm^2")
#    print(f"scattering angles: {SingleAngleRMSFromCum[:j+1]}")

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
sigmaRange = layerThickness/np.sqrt(2)*CumSigmaFit**2
plt.figure(figsize=(10, 7))
plotSingleThickness(20, depth, rangeDecrease, r"$z/\cos(\theta)-z$")
plt.plot(G4depths, sigmaRange, color="navy", linewidth=2, label="Geant4 Lateral Scattering RMS")
plt.plot(G4depths, rangeDecreaseSigma, color="navy", linewidth=2, label="rangeDecreaseSigma")
plt.xlabel(r"$\Delta X$ / degree")
plt.ylabel("Probability density")
plt.yscale("log")
plt.title(f"$\\Delta X$ distribution at depth = {20} mm")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()