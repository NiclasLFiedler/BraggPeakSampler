import numpy as np
import matplotlib.pyplot as plt
import json
import ROOT
from scipy.optimize import curve_fit

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

with open("../../analysis/config.json", "r") as file:
    fullConfig = json.load(file)

detectorSelect = fullConfig["detectorSelect"]
targetSelect   = fullConfig["targetSelect"]
plotEnable     = fullConfig["plotEnable"]
config         = fullConfig["detectors"][detectorSelect]

datasetSelect        = config["datasetSelect"]
detectorType         = config["detectorType"]
beamEnergy           = config["beamEnergy"]
nLayers              = config["nLayers"]
crystalSize          = config["crystalSize"]
gapSizeZ             = config["gapSizeZ"]
secondaryLayerStatus = config["secondaryLayerStatus"]
nSecondaryLayers     = config["nSecondaryLayers"]
secLayerSizeZ        = config["secLayerSizeZ"]
absorberStatus       = config["absorberStatus"]
absorberSize         = config["absorberSize"]
reversedStatus       = config["reversedStatus"]
normStatus           = config["normStatus"]
simulationStatus     = config["simulationStatus"]
teflonThickness      = config["teflonThickness"]
aluThickness         = config["aluThickness"]
coincidenceTime      = config["coincidenceTime"]
coincidenceLayer     = config["coincidenceLayer"]
targetThickness      = config["targetThickness"]

detectorZ = crystalSize[2]

in_data  = ["notarget", "homotarget", "heterotarget"]
datapath = f"data/{in_data[targetSelect]}"

file = ROOT.TFile("data/temp/data.root")
tree = file.Get("vtree")

# --- Accumulation dictionaries ---
energy_per_layer  = {}
energy_per_layer2 = {}
energy_per_event  = {}
wet_per_layer     = {}

for event in tree:
    eventid  = event.event
    layer    = event.NDet
    edep     = event.EDep
    wetAccum = event.WetAccum
    
    energy_per_layer[layer]  = energy_per_layer.get(layer, 0.0) + edep
    energy_per_layer2[layer] = energy_per_layer2.get(layer, 0.0) + edep
    energy_per_event[eventid] = energy_per_event.get(eventid, 0.0) + edep

    
    if layer not in wet_per_layer:
        wet_per_layer[layer] = []
    wet_per_layer[layer].append(wetAccum)

layers = np.array(sorted(energy_per_layer.keys()))
depth  = (layers + 0.5) * detectorZ / nLayers
energy = np.array([energy_per_layer[l]  for l in layers])
energy2 = np.array([energy_per_layer2[l] for l in layers])


depthWET = []
mean_wet_per_layer     = {}
mean_wet_err_per_layer = {}

for layer, wet_values in wet_per_layer.items():
    wet_values = np.array(wet_values)

    counts, bins = np.histogram(wet_values, bins=2000)
    bin_centers  = 0.5 * (bins[1:] + bins[:-1])

    max_idx      = np.argmax(counts[30:])
    max_idx += 30
    
    peak_pos     = bin_centers[max_idx] 
    fit_mask     = (bin_centers > peak_pos * 0.7) & (bin_centers < peak_pos * 1.3)
    x_fit, y_fit = bin_centers[fit_mask], counts[fit_mask]


    if len(x_fit) < 3:
        # Fallback to simple mean if not enough points to fit
        mean_wet_per_layer[layer]     = np.mean(wet_values)
        mean_wet_err_per_layer[layer] = np.std(wet_values) / np.sqrt(len(wet_values))
        continue

    try:
        p0   = [counts[max_idx], peak_pos, np.std(wet_values)]
        popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=p0, maxfev=1000000)
        mean_wet_per_layer[layer]     = popt[1]           # mu of Gaussian
        mean_wet_err_per_layer[layer] = np.sqrt(pcov[1,1])  # uncertainty on mu
    except RuntimeError:
        mean_wet_per_layer[layer]     = np.mean(wet_values)
        mean_wet_err_per_layer[layer] = np.std(wet_values) / np.sqrt(len(wet_values))

    if(layer == 320):
        plt.hist(wet_values, bins=1000, alpha=0.6, label=f"Layer {layer}")
        plt.plot(x_fit, gaussian(x_fit, *popt), color='r', label="Fit region")
        # plt.yscale("log")
        plt.xlabel("WET (mm)")
        plt.ylabel("Counts")
        plt.title(f"WET distribution for layer {layer}")
        # plt.ylim(1, counts.max() * 10)
        plt.legend()
        # plt.show()
        plt.close()
for index, value in enumerate(mean_wet_per_layer):
    if index == 0:
        depthWET.append(mean_wet_per_layer[value]/20)
    else:
        depthWET.append(depthWET[-1]+mean_wet_per_layer[value]/20)
# Convert to arrays sorted by layer
mean_wet     = np.array([mean_wet_per_layer[l]     for l in layers])
mean_wet_err = np.array([mean_wet_err_per_layer[l] for l in layers])

# --- Plot WET histograms per layer ---
n_cols = 4
n_rows = int(np.ceil(len(layers) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows))
axes = axes.flatten()

for idx, layer in enumerate(layers):
    wet_values  = np.array(wet_per_layer[layer])
    counts, bins = np.histogram(wet_values, bins=100)
    bin_centers  = 0.5 * (bins[1:] + bins[:-1])

    axes[idx].bar(bin_centers, counts, width=(bins[1]-bins[0]), alpha=0.6, label="WET")

    # Overlay Gaussian fit if successful
    if layer in mean_wet_per_layer:
        mu_fit  = mean_wet_per_layer[layer]
        fit_mask = (bin_centers > mu_fit * 0.9) & (bin_centers < mu_fit * 1.1)
        if fit_mask.sum() >= 3:
            x_plot = np.linspace(bin_centers[fit_mask][0], bin_centers[fit_mask][-1], 300)
            A_fit  = np.max(counts[fit_mask])
            sig_fit = np.std(wet_values)
            axes[idx].plot(x_plot, gaussian(x_plot, A_fit, mu_fit, sig_fit),
                           'r-', linewidth=1.5, label=f"μ={mu_fit:.2f} mm")

    axes[idx].set_title(f"Layer {layer}")
    axes[idx].set_xlabel("WET (mm)")
    axes[idx].set_ylabel("Counts")
    axes[idx].legend(fontsize=7)


plt.suptitle("WET distribution per layer", fontsize=14)
plt.tight_layout()
plt.savefig(f"{datapath}/input/wet_histograms_{nLayers}.png", dpi=150)
plt.close()

# --- Plot mean WET vs layer ---
plt.figure()
plt.errorbar(layers, mean_wet, yerr=mean_wet_err, fmt='o-', capsize=3)
plt.xlabel("Layer index")
plt.ylabel("Mean WET (mm)")
plt.title("Mean WET per layer")
plt.tight_layout()
plt.savefig(f"{datapath}/input/mean_wet_{nLayers}.png", dpi=150)
plt.close()

# --- Existing energy plots ---
plt.plot(depth, energy)
plt.xlabel("Depth (mm)")
plt.ylabel("Total Deposited Energy")
plt.title("Bragg Curve")
plt.close()

energy_values = np.array(list(energy_per_event.values()))
counts, bins  = np.histogram(energy_values, bins=1000)
bin_centers   = 0.5 * (bins[1:] + bins[:-1])
max_index     = np.argmax(counts[30:])
peak_position = bin_centers[max_index]
fit_mask      = bin_centers>10 & (bin_centers > peak_position * 0.9) & (bin_centers < peak_position * 1.1)
x_fit, y_fit  = bin_centers[fit_mask], counts[fit_mask]
print(f"Peak position: {peak_position}")
print(f"Fit range: {x_fit[0]} to {x_fit[-1]} with {len(x_fit)} points")
print(f"Counts in fit range: {y_fit.sum()} out of {counts.sum()} total counts")
print(f"Amplitude: {np.max(y_fit)}")
A0, mu0, sigma0 = np.max(y_fit), peak_position, peak_position * 0.02
popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[A0, mu0, sigma0], maxfev=1000000)
A, mu, sigma = popt
perr     = np.sqrt(np.diag(pcov))
mu_err   = perr[1]
sigma_err = perr[2]

print(f"Mu: {mu} +- {mu_err}")
print(f"Sigma: {sigma} +- {sigma_err}")

plt.hist(energy_values, bins=1000, alpha=0.6, label="Data")
x_plot = np.linspace(min(x_fit), max(x_fit), 500)
plt.plot(x_plot, gaussian(x_plot, *popt), 'r-', linewidth=2, label="Gaussian Fit")
plt.xlabel("Cumulative deposited energy")
plt.ylabel("Counts")
plt.legend()
plt.show()
# plt.close()

save_kwargs = dict(
    depth      = depth,
    depthWET  = depthWET,       
    depth_err  = (abs(depth[1] - depth[0])) / np.sqrt(12),
    dose       = energy,
    dose_err   = np.sqrt(energy),
    energy     = energy_values,
    amplitude  = A,
    mean       = mu,
    sigma      = sigma,
    mean_err   = mu_err,
    sigma_err  = sigma_err,
    mean_wet   = mean_wet,
    mean_wet_err = mean_wet_err,
)

if targetSelect == 0:
    np.savez(f"{datapath}/input/depthdose{nLayers}.npz", **save_kwargs)
else:
    np.savez(f"{datapath}/input/depthdose{nLayers}_{targetThickness}.npz", **save_kwargs)