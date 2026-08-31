import uproot
import numpy as np
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
    var = tree["CumVariance"].array(library="np")

theta = np.degrees(np.sqrt(var))

Vardepths = np.unique(depth)

mean_theta = np.array([
    np.mean(theta[depth == d])
    for d in Vardepths
])

E0 = 220  # Initial kinetic energy in MeV
m_p = 938.272  # Proton rest mass in MeV/c^2
X0 = 36.08   # Radiation length of water in cm

data = analysisFunctions.load_EnergyRange("../../range_energy/data_analysis/h2o_alt_range_energy.npz")

alpha = data.alpha[0]
p_exp = data.p

R0 = alpha * (E0**p_exp)

N = len(mean_theta)
Vardepths = Vardepths[:N]
mask = Vardepths < R0

Vardepths = Vardepths[mask]
mean_theta = mean_theta[mask]


print(f"R0: {R0}")
dx = 0.05
x_max = 0.999 * R0
depths = np.arange(0, x_max, dx)

E_k = E0 * np.maximum(0, (1 - depths / R0)) ** (1.0 / p_exp)
beta_p = (E_k * (E_k + 2 * m_p)) / (E_k + m_p)
integrand = (13.6 / beta_p) ** 2 * dx / X0
uncorrected_variance = np.cumsum(integrand)

theta_iH_rad = np.zeros_like(depths)

depthDiff = (depths-depths[0])
depthDiff[0] = depthDiff[1]
# print(depthDiff)
global_log_factor = 1 + 0.038 * np.log(depthDiff / X0)
# global_log_factor = 1 + 0.038 * np.log(depths[1:]/10 / X0)

theta_iH_rad[1:] = global_log_factor[1:] * np.sqrt(uncorrected_variance[1:])
theta_iH_deg = np.degrees(theta_iH_rad)


log_step = 1 + 0.038 * np.log(dx / X0)
theta_step = (13.6 / beta_p) * np.sqrt(dx / X0) * log_step
theta_naive_deg = np.degrees(np.sqrt(np.cumsum(theta_step**2)))

sigma_r = depths[-1]*(1-np.cos(theta_iH_deg[-1])) 
print(f"Depth of material: {depths[-1]:.4f} cm")
print(f"Final cumulative scattering angle: {theta_iH_deg[-1]:.4f} mrad")

print(f"Final range sigma: {sigma_r:.4f} cm")



plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.plot(depths, theta_iH_deg, color="navy", linewidth=2, label="Integral Highland (Thick Target)")
# plt.plot(depths, theta_naive_deg, color="black", linewidth=2.5, label="True Integral Highland (Thick Target)")
plt.scatter(
    Vardepths,
    mean_theta,
    marker="o",
    s=10,
    color="orange",
    label="Geant4"
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