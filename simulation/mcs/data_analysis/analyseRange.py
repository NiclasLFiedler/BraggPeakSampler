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
thetaX = np.degrees(np.arctan(1/thetaX))
thetaY = np.degrees(np.arctan(1/thetaY))
# --------------------------------------------------
# Define depth bins
# --------------------------------------------------

bin_width = 1.0  # mm

bins = np.arange(0, np.max(depth) + bin_width, bin_width)
centers = 0.5 * (bins[:-1] + bins[1:])

# --------------------------------------------------
# Collect one angle per event per depth bin
# --------------------------------------------------

anglesX = [[] for _ in range(len(centers))]
anglesY = [[] for _ in range(len(centers))]

unique_events = np.unique(event)

for evt in unique_events:

    mask = event == evt

    z_evt = depth[mask]
    tx_evt = thetaX[mask]
    ty_evt = thetaY[mask]

    # Sort steps by depth

    order = np.argsort(z_evt)
    
    z_evt = z_evt[order]
    tx_evt = tx_evt[order]
    ty_evt = ty_evt[order]

    # Which depth bin does each step belong to?
    bin_index = np.digitize(z_evt, bins) - 1

    # Keep only ONE value per depth bin
    seen = set()

    for i in range(len(z_evt)):

        b = bin_index[i]

        if b < 0 or b >= len(centers):
            continue

        if b not in seen:
            anglesX[b].append(tx_evt[i])
            anglesY[b].append(ty_evt[i])
            seen.add(b)

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
# Plot
# --------------------------------------------------

plt.figure(figsize=(8,5))

plt.plot(centers, sigmaX, label=r'$\sigma(\theta_x)$')
plt.plot(centers, sigmaY, label=r'$\sigma(\theta_y)$')

plt.xlabel("Depth (mm)")
plt.ylabel("RMS projected angle (mrad)")
plt.title("Multiple Coulomb Scattering")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()