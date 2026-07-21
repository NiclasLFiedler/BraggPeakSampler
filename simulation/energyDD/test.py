import matplotlib.pyplot as plt
import numpy as np

# ==========================================
# 1. Physical Constants & Parameters
# ==========================================
E0 = 220  # Initial kinetic energy in MeV
m_p = 938.272  # Proton rest mass in MeV/c^2
X0 = 36.08  # Radiation length of water in cm

# Empirical water range parameters (Bortfeld model)
alpha = 0.0022  # Range coefficient (cm / MeV^1.77)
p_exp = 1.77  # Range exponent

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

# ==========================================
# 6. Plotting Results
# ==========================================
fig, ax = plt.subplots(figsize=(9, 5.5))

# Plot True Integral Highland

sigma_r = depths[-1]/10*(1-np.cos(theta_iH_deg[-1])) 
print(f"Depth of material: {depths[-1]:.4f} cm")
print(f"Final cumulative scattering angle: {theta_iH_deg[-1]:.4f} mrad")

print(f"Final range sigma: {sigma_r:.4f} cm")

ax.plot(
    depths,
    theta_iH_deg,
    color="navy",
    linewidth=2.5,
    label="True Integral Highland (Thick Target)",
)

# Plot Naive Local Highland for reference
ax.plot(
    depths,
    theta_naive_deg,
    color="darkred",
    linestyle="--",
    linewidth=1.8,
    label="Naive Step-Wise Highland (Local Log Term)",
)

# Vertical line highlighting the stopping range
ax.axvline(
    x=R0,
    color="gray",
    linestyle=":",
    linewidth=1.5,
    label=f"Proton Range $R_0 \\approx {R0:.2f}$ cm",
)

ax.set_xlabel("Depth in Water (cm)", fontsize=12)
ax.set_ylabel("Cumulative RMS Angle $\Theta_{\mathrm{cum}}$ (degrees)", fontsize=12)
ax.set_title(
    f"Proton Multiple Coulomb Scattering in Water\n"
    f"Integral Highland (Thick Target) vs. Naive Formula ($E_0 = {E0}$ MeV)",
    fontsize=13,
    pad=15,
)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(fontsize=11)

fig.tight_layout()
plt.show()