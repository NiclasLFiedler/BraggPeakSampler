import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Define the exponential decay function ---
def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

# --- Load histogram from ROOT file ---
file_path = "c_sipm_3.root"
hist_name = "h_ampl0"

with uproot.open(file_path) as file:
    hist = file[hist_name]
    bin_edges = hist.axis().edges()
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    counts = hist.values()

# --- Define fit range (tail only, e.g. last 20% of nonzero bins) ---
#nonzero_indices = np.nonzero(counts)[0]
#tail_start = int(nonzero_indices[-1] - 0.5 * len(nonzero_indices))
tail_start = int(len(bin_centers) * 0.45)
nonzero_indices = np.nonzero(counts)[0]
print("Nonzero indices:", len(nonzero_indices))

tail_start = int(nonzero_indices[-1] - 0.8 * len(nonzero_indices))

print(len(counts))
print(tail_start) 
x_tail = bin_centers[tail_start:]
y_tail = counts[tail_start:]

# --- Fit the exponential decay to the tail ---
popt, pcov = curve_fit(exp_decay, x_tail, y_tail, p0=(1e7, 2000), maxfev=1000000000)

mask = exp_decay(x_tail, *popt) > 0.1

# --- Plot the histogram with logarithmic y-axis ---
plt.figure(figsize=(10, 6))
plt.plot(bin_centers, counts, drawstyle='steps-mid', label="Histogram")
plt.plot(x_tail[mask], exp_decay(x_tail[mask], *popt), 'r--', label=f"Fit: A*exp(-x/τ)\nA={popt[0]:.2e}, τ={popt[1]:.2f}")
plt.yscale("log")
plt.xlabel("x")
plt.ylabel("Counts")
plt.title("Histogram with Exponential Fit")
plt.legend()
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()