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
    """
    Calculate Gaussian-fit sigma for every requested depth.
    """

    sigma_fit = np.full(len(depths), np.nan)
    variance_fit = np.full(len(depths), np.nan)

    std_data = np.full(len(depths), np.nan)
    variance_data = np.full(len(depths), np.nan)

    for i, d in enumerate(depths):

        selected = angles[
            np.abs(depth - d) < tolerance
        ]

        selected = selected[np.isfinite(selected)]

        if len(selected) < 10:
            continue

        counts, bin_edges = np.histogram(
            selected,
            bins=bins,
            density=True
        )

        bin_centers = (
            0.5 * (bin_edges[:-1] + bin_edges[1:])
        )

        A0 = np.max(counts)
        mu0 = np.mean(selected)
        sigma0 = np.std(selected)

        try:
            popt, pcov = curve_fit(
                gaussian,
                bin_centers,
                counts,
                p0=[A0, mu0, sigma0]
            )

            A, mu, sigma = popt

            sigma = abs(sigma)

            sigma_fit[i] = sigma
            variance_fit[i] = sigma**2

        except RuntimeError:
            continue

        std_data[i] = np.std(selected)
        variance_data[i] = np.var(selected)

    return (
        sigma_fit,
        variance_fit,
        std_data,
        variance_data
    )

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

##########test
CumSigmaFit, _, CumStd, _ =  gaussian_sigma_vs_depth(
        depth,
        CumScatteringAngle,
        G4depths,
        tolerance=0.1,
        bins=2000
    )

plt.figure(figsize=(12, 8))

plt.plot(
    G4depths,
    CumSigmaFit,
    "o-",
    label="Gaussian fit σ"
)

plt.plot(
    G4depths,
    CumAngleRMS,
    "s--",
    label=r"$\sqrt{\mathrm{np.var}}$"
)

plt.xlabel("Depth / mm")
plt.ylabel("Cumulative scattering angle RMS / degree")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
##########test end

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
    N = len(CumAngleRMS)
    G4depths = G4depths[:N]
    SingleAngleRMS = SingleAngleRMS[:N]

    mask = G4depths < R0 + 1

    G4depths = G4depths[mask]
    CumAngleRMS = CumAngleRMS[mask]
    SingleAngleRMS = SingleAngleRMS[mask]

layerThickness  =   1
lateralVariance = np.zeros(len(G4depths))

for j, d in enumerate(G4depths):

    lever_arm = d - G4depths[:j] - layerThickness / 2

    lateralVariance[j] = np.sum((lever_arm * np.radians(SingleAngleRMS[:j]))**2)

m_p = 938.272
X0 = 36.08
dx = 0.01
x_max = 0.999 * R0
depths = np.arange(dx, x_max, dx)

E_k = E0 * (1 - depths/R0)**(1/p_exp)
beta_p = E_k * (E_k + 2*m_p)/(E_k + m_p)

integrand = (13.6/beta_p)**2 * dx/X0
uncorrected_variance = np.cumsum(integrand)
log_factor = 1 + 0.038*np.log(depths/X0)
theta_integrated_rad = log_factor * np.sqrt(uncorrected_variance)
theta_integrated_deg = np.degrees(theta_integrated_rad)

# plt.figure(figsize=(12, 9))
# plt.plot(depths, theta_integrated_rad*1000, linewidth=2, label=r"$\beta p$")

# plt.xlabel("Depth / cm")
# plt.ylabel(r"$\beta p$ / MeV")
# plt.title(r"Proton $\beta p$ as a function of depth")
# plt.grid(True)
# plt.legend()

# plt.tight_layout()
# plt.show()

log_step = 1 + 0.038 * np.log(dx / X0)
theta_step = (13.6 / beta_p) * np.sqrt(dx / X0) * log_step
theta_naive_deg = np.degrees(np.sqrt(np.cumsum(theta_step**2)))

print(f"Depth of material: {depths[-1]:.4f} cm")
print(f"Final cumulative scattering angle: {theta_integrated_deg[-1]:.4f} mrad")

plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(12, 9))

plt.plot(depths, theta_integrated_deg, color="navy", linewidth=2, label="Integral Highland (Thick Target)")
plt.plot(depths, theta_naive_deg, color="black", linewidth=2.5, label="Naive Integral Highland (Thick Target)")
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
    CumSigmaFit,
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