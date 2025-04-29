import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Define the exponential decay function ---
def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

# --- Load histogram from ROOT file ---
file_path = "c_sipm_1.root"
hist_name = "h_ampl0"

with uproot.open(file_path) as file:
    hist = file[hist_name]
    bin_edges = hist.axis().edges()
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    counts = hist.values()

# --- Define fit range (tail only, e.g. last 20% of nonzero bins) ---
nonzero_indices = np.nonzero(counts)[0]
tail_start = int(nonzero_indices[-1] - 0.5 * len(nonzero_indices))
x_tail = bin_centers[tail_start:]
y_tail = counts[tail_start:]

# --- Fit the exponential decay to the tail ---
popt, pcov = curve_fit(exp_decay, x_tail, y_tail, p0=(y_tail[0], 10))

# --- Plot the histogram with logarithmic y-axis ---
plt.figure(figsize=(10, 6))
plt.plot(bin_centers, counts, drawstyle='steps-mid', label="Histogram")
plt.plot(x_tail, exp_decay(x_tail, *popt), 'r--', label=f"Fit: A*exp(-x/τ)\nA={popt[0]:.2e}, τ={popt[1]:.2f}")
plt.yscale("log")
plt.xlabel("x")
plt.ylabel("Counts")
plt.title("Histogram with Exponential Fit")
plt.legend()
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()