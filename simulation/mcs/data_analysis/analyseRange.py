import uproot
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# Load ROOT file
# --------------------------------------------------

file = uproot.open("h2oproj.root")
tree = file["braggsampler"]

event  = tree["event"].array(library="np")
depth  = tree["pos"].array(library="np")      # z-position in mm
thetaX = tree["thetaX"].array(library="np")
thetaY = tree["thetaY"].array(library="np")
thetaX = np.degrees(np.arctan(thetaX))
thetaY = np.degrees(np.arctan(thetaY))

# --------------------------------------------------
# Define depth bins
# --------------------------------------------------

bin_width = 0.5  # mm
bins = np.arange(0, np.max(depth) + bin_width, bin_width)
centers = 0.5 * (bins[:-1] + bins[1:])

# --------------------------------------------------
# Vectorized: Assign depth bins and sort by event+depth
# --------------------------------------------------

# Assign each measurement to a depth bin
bin_index = np.digitize(depth, bins) - 1

# Filter out out-of-bounds
valid_mask = (bin_index >= 0) & (bin_index < len(centers))
event = event[valid_mask]
depth = depth[valid_mask]
thetaX = thetaX[valid_mask]
thetaY = thetaY[valid_mask]
bin_index = bin_index[valid_mask]

# Sort by event ID, then by depth (to get sorted steps within each event)
sort_order = np.lexsort((depth, event))
event = event[sort_order]
depth = depth[sort_order]
thetaX = thetaX[sort_order]
thetaY = thetaY[sort_order]
bin_index = bin_index[sort_order]

# --------------------------------------------------
# Find event boundaries and keep ONE value per event per bin
# --------------------------------------------------

# Identify where event changes
event_changes = np.concatenate(([True], event[1:] != event[:-1], [True]))
event_boundaries = np.where(event_changes)[0]

# Initialize output arrays for all bins
anglesX = [[] for _ in range(len(centers))]
anglesY = [[] for _ in range(len(centers))]

# Process each event
for i in range(len(event_boundaries) - 1):
    start_idx = event_boundaries[i]
    end_idx = event_boundaries[i + 1]
    
    # Get data for this event
    event_bins = bin_index[start_idx:end_idx]
    event_thetaX = thetaX[start_idx:end_idx]
    event_thetaY = thetaY[start_idx:end_idx]
    
    # Get unique bins (preserves first occurrence due to sorting)
    unique_bins, first_indices = np.unique(event_bins, return_index=True)
    
    # Append first occurrence of each bin
    for bin_id, idx in zip(unique_bins, first_indices):
        anglesX[bin_id].append(event_thetaX[idx])
        anglesY[bin_id].append(event_thetaY[idx])

# --------------------------------------------------
# Calculate RMS
# --------------------------------------------------

sigmaX = np.array([
    np.std(a, ddof=1) if len(a) > 1 else np.nan
    for a in anglesX
])

sigmaY = np.array([
    np.std(a, ddof=1) if len(a) > 1 else np.nan
    for a in anglesY
])

# --------------------------------------------------
# Plot 1: RMS vs Depth
# --------------------------------------------------

# ==========================================
# 1. Physical Constants & Parameters
# ==========================================
E0 = 210  # Initial kinetic energy in MeV
m_p = 938.272  # Proton rest mass in MeV/c^2
X0 = 36.08  # Radiation length of water in cm

# Empirical water range parameters (Bortfeld model)
alpha = 0.02369  # Range coefficient (cm / MeV^1.77)
p_exp = 1.757  # Range exponent

# Total projected range in water (cm)
R0 = alpha * (E0**p_exp)

# ==========================================
# 2. Simulation Grid Setup
# ==========================================
dx = 0.005  # Step size in cm (0.05 mm)
x_max = 0.999 * R0  # Stop just short of R0 to avoid division by zero
depths = np.arange(0, x_max, dx)

# ==========================================
# 3. Energy Loss & Kinematics Along Track
# ==========================================
# Kinetic energy E_k(x) as a function of depth
E_k = E0 * np.maximum(0, (1 - depths / R0)) ** (1.0 / p_exp)

# Kinematic factor beta * c * p (in MeV)
beta_p = (E_k * (E_k + 2 * m_p)) / (E_k + m_p)

# ==========================================
# 4. True Integral Highland (Thick Target)
# ==========================================
# Uncorrected scattering integrand f(x') = (13.6 / beta_p)^2 / X0
integrand = (13.6 / beta_p) ** 2 / X0

# Numerical integration along the path: \int_0^x f(x') dx'
uncorrected_variance = np.cumsum(integrand * dx)

# Initialize output array
theta_iH_rad = np.zeros_like(depths)

# Global logarithmic correction factor for depth x > 0: [1 + 0.038 * ln(x / X0)]
# (Skipping x=0 to prevent ln(0) errors)
global_log_factor = 1 + 0.038 * np.log(depths[1:] / X0)

# True Integral Highland calculation
theta_iH_rad[1:] = global_log_factor * np.sqrt(uncorrected_variance[1:])
theta_iH_deg = np.degrees(theta_iH_rad)

# ==========================================
# 5. Naive Step-Wise Highland (For Comparison)
# ==========================================
# Applies the log factor locally per dx step (step-size dependent)
log_step = 1 + 0.038 * np.log(dx / X0)
theta_step = (13.6 / beta_p) * np.sqrt(dx / X0) * log_step
theta_naive_deg = np.degrees(np.sqrt(np.cumsum(theta_step**2)))

sigma_r = depths[-1]/10*(1-np.cos(theta_iH_deg[-1])) 
print(f"Depth of material: {depths[-1]:.4f} cm")
print(f"Final cumulative scattering angle: {theta_iH_deg[-1]:.4f} mrad")

print(f"Final range sigma: {sigma_r:.4f} cm")



plt.figure(figsize=(8, 5))

plt.plot(depths, theta_iH_deg, color="navy", linewidth=2.5, label="True Integral Highland (Thick Target)")
plt.plot(centers, sigmaX, label=r'$\sigma(\theta_x)$')
plt.plot(centers, sigmaY, label=r'$\sigma(\theta_y)$')

plt.xlabel("Depth (mm)")
plt.ylabel("RMS projected angle (mrad)")
plt.title("Multiple Coulomb Scattering")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

# --------------------------------------------------
# Plot 2: Angle distributions at different depths
# --------------------------------------------------

n_plots = 10
indices = np.linspace(0, len(centers) - 1, n_plots, dtype=int)

plt.figure(figsize=(10, 8))

for idx in indices:
    angles = np.array(anglesX[idx])

    if len(angles) < 10:
        continue

    plt.hist(
        angles,
        bins=500,
        histtype="step",
        density=True,
        label=f"{centers[idx]:.1f} mm"
    )

plt.xlabel(r"$\theta_x$ (mrad)")
plt.ylabel("Normalized counts")
plt.yscale("log")
plt.title("Projected MCS angle distributions at different depths")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()