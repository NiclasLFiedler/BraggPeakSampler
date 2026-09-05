import uproot
import numpy as np
import matplotlib
matplotlib.use("QtAgg")
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
# --------------------------------------------------
# Load ROOT file
# --------------------------------------------------

with uproot.open("h2oproj.root") as f:
    tree = f["braggsampler"]

    event = tree["event"].array(library="np")
    depth = tree["depth"].array(library="np")
    CumScatteringAngle =   np.degrees(tree["CumScatteringAngle"].array(library="np"))
    ScatteringAngle =   np.degrees(tree["SingleScatteringAngle"].array(library="np"))

G4depths = np.unique(depth)

CumVariance = np.array([np.var(CumScatteringAngle[depth == d]) for d in G4depths])
ScatteringVariance= np.array([np.var(ScatteringAngle[depth == d]) for d in G4depths])

CumAngleRMS = np.sqrt(CumVariance)
SingleAngleRMS = np.sqrt(ScatteringVariance)

plt.figure(figsize=(12, 8))

plt.plot( G4depths, CumAngleRMS, "s--", markersize=4, label="Cumulative RMS")
plt.plot( G4depths, SingleAngleRMS, "^-", markersize=4, label="Single scattering angle RMS")

plt.xlabel("Depth / mm")
plt.ylabel("Scattering angle RMS / degree")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.close()
# plt.show()

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

useMask = True
if useMask:
    N = len(CumAngleRMS)
    G4depths = G4depths[:N]
    SingleAngleRMS = SingleAngleRMS[:N]

    mask = G4depths < R0 + 1

    G4depths = G4depths[mask]
    CumAngleRMS = CumAngleRMS[mask]
    SingleAngleRMS = SingleAngleRMS[mask]

layerThickness  = 0.1
lateralVariance = [(layerThickness/2*np.deg2rad(angle))**2 for index, angle in enumerate(SingleAngleRMS)]
cumulativeVariance = np.cumsum(lateralVariance)

m_p = 938.272
X0 = 36.08
dx = 0.1
x_max = 0.999 * R0
depths = np.arange(0, x_max, dx)

E_k = E0 * np.maximum(0, (1 - depths / R0)) ** (1.0 / p_exp)
beta_p = (E_k * (E_k + 2 * m_p)) / (E_k + m_p)
integrand = (13.6 / beta_p) ** 2 * dx / X0
uncorrected_variance = np.cumsum(integrand)

theta_iH_rad = np.zeros_like(depths)

depthDiff = (depths-depths[0])
depthDiff[0] = depthDiff[1]

global_log_factor = 1 + 0.038 * np.log(depthDiff / X0)
# global_log_factor = 1 + 0.038 * np.log(depths[1:]/10 / X0)

theta_iH_rad[1:] = global_log_factor[1:] * np.sqrt(uncorrected_variance[1:])
theta_iH_deg = np.degrees(theta_iH_rad)


log_step = 1 + 0.038 * np.log(dx / X0)
theta_step = (13.6 / beta_p) * np.sqrt(dx / X0) * log_step
theta_naive_deg = np.degrees(np.sqrt(np.cumsum(theta_step**2)))

E_k = E0 * np.maximum(0, (1 - depths / R0)) ** (1.0 / p_exp)

beta_p = np.sqrt(E_k * (E_k + 2 * m_p)) / (E_k + m_p)

t_over_X0 = depths / X0

theta_highland_rad = (13.6 / beta_p) * np.sqrt(t_over_X0) * (1 + 0.038 * np.log(t_over_X0))

theta_highland_deg = np.degrees(theta_highland_rad)

variance_highland = theta_highland_rad ** 2


print(f"Depth of material: {depths[-1]:.4f} cm")
print(f"Final cumulative scattering angle: {theta_iH_deg[-1]:.4f} mrad")



plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.plot(depths, theta_highland_deg, color="navy", linewidth=2, label="Integral Highland (Thick Target)")
# plt.plot(depths, theta_naive_deg, color="black", linewidth=2.5, label="True Integral Highland (Thick Target)")
plt.scatter(
    G4depths,
    CumAngleRMS,
    marker="o",
    s=10,
    color="orange",
    label="Geant4"
)

plt.scatter(
    G4depths,
    SingleAngleRMS,
    marker="o",
    s=10,
    color="green",
    label="Geant4 Theta"
)


plt.xlabel("Depth / cm")
plt.ylabel("RMS projected angle / degree")
plt.title("Multiple Coulomb Scattering")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig(
    "multiple_coulomb_scattering.svg",
    format="svg",
    bbox_inches="tight"
)
plt.show()


plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.scatter(
    G4depths,
    lateralVariance,
    marker="o",
    s=10,
    color="orange",
    label="Geant4"
)

plt.scatter(
    G4depths,
    np.sqrt(cumulativeVariance),
    marker="o",
    s=10,
    color="green",
    label="Geant4"
)

plt.xlabel("Depth / cm")
plt.ylabel("Lateral Scattering / cm")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig(
    "lateral_scattering.svg",
    format="svg",
    bbox_inches="tight"
)
plt.show()