import numpy as np
import matplotlib.pyplot as plt
import json
import ROOT

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

for event in tree:
    layer = event.NDet
    edep  = event.EDep
    
    if layer not in energy_per_layer:
        energy_per_layer[layer] = 0.0
        
    energy_per_layer[layer] += edep

# Convert to arrays
layers = np.array(sorted(energy_per_layer.keys()))
depth = (layers+0.5)*50/nLayers
energy = np.array([energy_per_layer[l] for l in layers])

plt.plot(depth, energy)
plt.xlabel("Detector Layer (NDet)")
plt.ylabel("Total Deposited Energy")
plt.title("Bragg Curve")
plt.show()

if(targetSelect == 0):
    np.savez(
        f"{datapath}/input/depthdose{nLayers}.npz",
        depth = depth,
        depth_err = (abs(depth[1]-depth[0]))/np.sqrt(12),
        dose = energy,
        dose_err = np.sqrt(energy)
        )
else:
    np.savez(
        f"{datapath}/input/depthdose{nLayers}_{targetThickness}.npz",
        depth = depth,
        depth_err = (abs(depth[1]-depth[0]))/np.sqrt(12),
        dose = energy,
        dose_err = np.sqrt(energy)
        )