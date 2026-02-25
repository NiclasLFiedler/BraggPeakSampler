import numpy as np
import matplotlib.pyplot as plt
import json
import ROOT
from scipy.optimize import curve_fit
import numpy as np

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

with open("../../analysis/config.json", "r") as file:
    fullConfig = json.load(file)

detectorSelect = fullConfig["detectorSelect"]
targetSelect = fullConfig["targetSelect"]
plotEnable = fullConfig["plotEnable"]

config = fullConfig["detectors"][detectorSelect]

datasetSelect = config["datasetSelect"]
detectorType = config["detectorType"]
beamEnergy = config["beamEnergy"]
nLayers = config["nLayers"]
crystalSize = config["crystalSize"]
gapSizeZ = config["gapSizeZ"]
secondaryLayerStatus = config["secondaryLayerStatus"]
nSecondaryLayers = config["nSecondaryLayers"]
secLayerSizeZ = config["secLayerSizeZ"]
absorberStatus = config["absorberStatus"]
absorberSize = config["absorberSize"]
reversedStatus = config["reversedStatus"]
normStatus = config["normStatus"]
simulationStatus = config["simulationStatus"]
teflonThickness = config["teflonThickness"]
aluThickness = config["aluThickness"]
coincidenceTime = config["coincidenceTime"]
coincidenceLayer = config["coincidenceLayer"]
targetThickness = config["targetThickness"]

in_data = ["notarget", "homotarget", "heterotarget"]
datapath = f"data/{in_data[targetSelect]}"

file = ROOT.TFile("data/temp/data.root")
tree = file.Get("vtree")

# Dictionary to accumulate energy per layer
energy_per_layer = {}
energy_per_layer2 = {}
energy_per_event = {}

for event in tree:
    eventid  = event.event
    layer = event.NDet
    edep  = event.EDep
    
    if layer not in energy_per_layer:
        energy_per_layer[layer] = 0.0
        
    energy_per_layer[layer] += edep
    
    energy_per_layer2[layer] = energy_per_layer2.get(layer, 0.0) + edep

    energy_per_event[eventid] = energy_per_event.get(eventid, 0.0) + edep

layers = np.array(sorted(energy_per_layer.keys()))
depth = (layers+0.5)*50/nLayers
energy = np.array([energy_per_layer[l] for l in layers])
energy2 = np.array([energy_per_layer2[l] for l in layers])

plt.plot(depth, energy)
plt.xlabel("Detector Layer (NDet)")
plt.ylabel("Total Deposited Energy")
plt.title("Bragg Curve")
plt.show()


energy_values = np.array(list(energy_per_event.values()))
counts, bins = np.histogram(energy_values, bins=2000)
bin_centers = 0.5 * (bins[1:] + bins[:-1])
max_index = np.argmax(counts)
peak_position = bin_centers[max_index]

fit_mask = (bin_centers > peak_position * 0.9) & \
           (bin_centers < peak_position * 1.1)

x_fit = bin_centers[fit_mask]
y_fit = counts[fit_mask]
# Initial guesses
A0 = np.max(y_fit)
mu0 = peak_position
sigma0 = peak_position * 0.02   # 2% guess

popt, pcov = curve_fit(gaussian, x_fit, y_fit,
                       p0=[A0, mu0, sigma0])

A, mu, sigma = popt
perr = np.sqrt(np.diag(pcov))

mu_err = perr[1]
sigma_err = perr[2]

print(f"Mu: {mu} +- {mu_err}")
print(f"Sigma: {sigma} +- {sigma_err}")

plt.hist(energy_values, bins=2000, alpha=0.6, label="Data")

x_plot = np.linspace(min(x_fit), max(x_fit), 500)
plt.plot(x_plot, gaussian(x_plot, *popt),
         'r-', linewidth=2, label="Gaussian Fit")

plt.xlabel("Cumulative deposited energy")
plt.ylabel("Counts")
plt.legend()
plt.show()

if(targetSelect == 0):
    np.savez(
        f"{datapath}/input/depthdose{nLayers}.npz",
        depth = depth,
        depth_err = (abs(depth[1]-depth[0]))/np.sqrt(12),
        dose = energy,
        dose_err = np.sqrt(energy),
        energy = energy_values,
        amplitude = A,
        mean = mu,
        sigma = sigma,
        mean_err = mu_err,
        sigma_err = sigma_err
        )
else:
    np.savez(
        f"{datapath}/input/depthdose{nLayers}_{targetThickness}.npz",
        depth = depth,
        depth_err = (abs(depth[1]-depth[0]))/np.sqrt(12),
        dose = energy,
        dose_err = np.sqrt(energy),
        energy = energy_values,
        amplitude = A,
        mean = mu,
        sigma = sigma,
        mean_err = mu_err,
        sigma_err = sigma_err
        )