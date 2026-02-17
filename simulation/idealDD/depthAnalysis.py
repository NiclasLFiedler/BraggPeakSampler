import numpy as np
import matplotlib.pyplot as plt
import json

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

meshSize = 50
sizeZ1 = meshSize
# sizeZ2 = 18
# sizeZ3 = meshSize-sizeZ1-sizeZ2

SegNumber1 = nLayers

# SegNumber2 = 800
# SegNumber3 = 300

bin_width1 = sizeZ1/SegNumber1
# bin_width2 = sizeZ2 / SegNumber2
# bin_width3 = sizeZ3 / SegNumber3

data1 = np.loadtxt(
    "build/meshEdep1.txt",
    delimiter=",",
    comments="#"
)


# data2 = np.loadtxt(
#     "build/meshEdep2.txt",
#     delimiter=",",
#     comments="#"
# )

# data3 = np.loadtxt(
#     "build/meshEdep3.txt",
#     delimiter=",",
#     comments="#"
# )

iZ1 = data1[:, 2]
edep1 = data1[:, 3]

# iZ2 = data2[:, 2]
# edep2 = data2[:, 3]*1/bin_width2

# iZ3 = data3[:, 2]
# edep3 = data3[:, 3]*1/bin_width3


depth1 = (abs(iZ1-SegNumber1)-0.5) * bin_width1
# depth2 = depth1[0] + (abs(iZ2-SegNumber2-1) + 0.5) * bin_width2
# depth3 = depth2[0] + (abs(iZ3-SegNumber3-1) + 0.5) * bin_width3
print(len(abs(iZ1-SegNumber1)-0.5))

# ---- Plot ----
plt.figure()
plt.step(depth1, edep1)
# plt.step(depth2, edep2)
# plt.step(depth3, edep3)

plt.xlabel("Depth in water (cm)")
plt.ylabel("Deposited energy (MeV)")
plt.title("Proton Depth Dose Curve")
plt.grid()

plt.close()

np.savez(
    f"{datapath}/input/depthdose{nLayers}.npz",
    depth = depth1,
    depth_err = (abs(depth1[-1]-depth1[-2]))/np.sqrt(12),
    dose = edep1,
    dose_err = np.sqrt(edep1)
)