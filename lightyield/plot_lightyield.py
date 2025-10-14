import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import ROOT

def gaussian(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

datapath = "data/1310_pedastal"
input_file = f"{datapath}/bps_spe.his.txt"
output_file = f"{datapath}/bps_spe.his.root"
hist_name = "hData"
hist_title = "Histogram from text file"
x_title = "x"
y_title = "Counts"

# --- Read data from text file ---
x_vals, y_vals = [], []
with open(input_file, "r") as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        x, y = map(float, parts[:2])
        x_vals.append(x)
        y_vals.append(y)

# --- Create histogram ---
nbins = len(x_vals)
x_min = min(x_vals)
x_max = max(x_vals)
hist = ROOT.TH1D(hist_name, f"{hist_title};{x_title};{y_title}", nbins, x_min, x_max)

# --- Fill histogram ---
for x, y in zip(x_vals, y_vals):
    bin_idx = hist.FindBin(x)
    hist.SetBinContent(bin_idx, y)

# --- Save histogram to ROOT file ---
root_file = ROOT.TFile(output_file, "RECREATE")
hist.Write()
root_file.Close()

print(f"\n✅ Histogram saved as '{hist_name}' in '{output_file}'.")