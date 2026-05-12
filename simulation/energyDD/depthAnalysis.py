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

for layer, wet_values in wet_per_layer.items():
    wet_values = np.array(wet_values)

    counts, bins = np.histogram(wet_values, bins=2000)
    bin_centers  = 0.5 * (bins[1:] + bins[:-1])

    max_idx      = np.argmax(counts)
    peak_pos     = bin_centers[max_idx]
    
    if(peak_pos <= 0.1):
        mean_wet_per_layer[layer] = mean_wet_per_layer[0]
    else:
        mean_wet_per_layer[layer] = peak_pos



    if(layer == 20):
        plt.hist(wet_values, bins=1000, alpha=0.6, label=f"Layer {layer}")
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
        depthWET.append(depthWET[-1]+(mean_wet_per_layer[value-1]+mean_wet_per_layer[value])/20)
# Convert to arrays sorted by layer
mean_wet     = np.array([mean_wet_per_layer[l] for l in layers])

# --- Plot mean WET vs layer ---
plt.figure()
plt.plot(layers, mean_wet)
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
x_plot = np.linspace(min(x_fit), max(x_fit), 2000)
plt.plot(x_plot, gaussian(x_plot, *popt), 'r-', linewidth=2, label="Gaussian Fit")
plt.xlabel("Cumulative deposited energy")
plt.ylabel("Counts")
plt.legend()
# plt.show()
plt.close()

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
)

if targetSelect == 0:
    np.savez(f"{datapath}/input/depthdose{nLayers}.npz", **save_kwargs)
else:
    np.savez(f"{datapath}/input/depthdose{nLayers}_{targetThickness}.npz", **save_kwargs)