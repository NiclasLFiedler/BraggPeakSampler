import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scienceplots
import uproot
from scipy.stats import truncnorm
import time
import pandas as pd
import ROOT
from scipy import special
import math
from scipy.signal import convolve
from scipy.signal import fftconvolve
from scipy.stats import norm
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy import ndimage
from scipy.integrate import simpson
import json

#plt.style.use(['science','notebook','grid']) 

a_h2o=2.585e-3
p_h2o=1.738
a_pwo=7.275e-4
p_pwo=1.690
a_dsb = 1.030e-3
p_dsb = 1.713

class fit_params:
    def __init__(self, Phi0=0, R0=0, sigma=0, epsilon=0, curve=[], stddev = []) -> None:
        self.Phi0 = Phi0
        self.R0 = R0
        self.sigma = sigma
        self.epsilon = epsilon
        self.curve = curve
        self.stddev = stddev

class fit_params_conv:
    def __init__(self, t=0, sigma=0, curve=[], cov = []) -> None:
        self.t = t
        self.sigma = sigma
        self.curve = curve
        self.cov = cov
def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-(x-mean)**2 / (2 * stddev**2)) / (np.sqrt(2 * np.pi) * stddev)

def right_sided_convolution(f, g, z_values):
    def convole(z):
        if(len(z_values)<100):
            zlin = np.linspace(0, 40, 1000)
            combined = np.concatenate([zlin, z_values])
            zFull = np.sort(combined)
        else:
            zFull = z_values

        indices = [np.where(zFull == value)[0][0] for value in z_values]

        dz = zFull[1]-zFull[0]
        convolution_result = np.zeros_like(zFull)
        
        for i, z_val in enumerate(zFull):
            shifted_f = lambda z: f(z + z_val)
            product = g(zFull) * shifted_f(zFull)
            convolution_result[i] = simpson(y=product, x=zFull)

        return convolution_result[indices]

    return convole(z_values)

def gaussian_with_cutoff(mean, sigma, cutoff=2.5):
    while True:
        value = np.random.normal(mean, sigma)
        if mean-cutoff <= value <= mean+cutoff:
            return value

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

datasets = ["MIT_05_2024", "simulation", "data"]
in_data = ["notarget", "homotarget", "heterotarget"]
in_title = ["without a target", "with the homogeneous target", "with the heterogeneous target"]

lineWidth = 2
capSize = 3

dataset = datasets[datasetSelect]
file = in_data[0]
targetFile = in_data[targetSelect]
if(targetSelect == 1):
    targetFile = f"{targetFile}{targetThickness}"

nosave = False

bhetero = False
if targetSelect == 2: bhetero = True

nbOfFits = 0
nbOfFitsHetero = 0
lineWidth = 2
capSize = 3

depth = []
depthTarget = []
depthErr = []
depthTargetErr = []
Dose = []
DoseTarget = []
DoseErr = []
DoseTargetErr = []
unfoldedEntries = []
unfoldedEntriesTarget = []


dataFile = np.load(f"{dataset}/{file}/input/depthdose{nLayers}.npz")
print(f"File: {dataFile}")

depth = np.flip(dataFile["depth"])
depthErr = np.flip(dataFile["depth_err"]) 
Dose = np.flip(dataFile["dose"]) 
DoseErr = np.flip(dataFile["dose_err"]) 
unfoldedEntries = dataFile["dose"][0]

if(bhetero):
    dataFile = np.load(f"{dataset}/{targetFile}/input/depthdose{nLayers}_{targetThickness}.npz")
    print(f"Hetero File: {dataFile}")
    depthTarget = np.flip(dataFile["depth"]) 
    depthTargetErr = np.flip(dataFile["depth_err"]) 
    DoseTarget = np.flip(dataFile["dose"]) 
    DoseTargetErr = np.flip(dataFile["dose_err"]) 
    unfoldedEntriesTarget = dataFile["dose"][0]

targetThicknesses = [52, 53, 54, 55, 56, 57]
if(targetSelect == 1):
    for thickness in targetThicknesses:
        dataFile = np.load(f"{dataset}/{targetFile}/input/depthdose{nLayers}_{targetThickness}.npz")
        depthTarget.append(np.flip(dataFile["depth"]))
        depthTargetErr.append(np.flip(dataFile["depth_err"]))
        DoseTarget.append(np.flip(dataFile["dose"]))
        DoseTargetErr.append(np.flip(dataFile["dose_err"]))
        unfoldedEntriesTarget.append(dataFile["dose"][0])

beta = 0.012
R0 = 30.99
if(targetSelect == 2):
    doseConversion = unfoldedEntries/unfoldedEntriesTarget*(1+beta*(R0-targetThickness))/(1+beta*R0)
    doseConversion = 1
    print(doseConversion)
    
    DoseTarget = np.array(DoseTarget)
    DoseTargetErr = np.array(DoseTargetErr)
    DoseTarget = DoseTarget * doseConversion
    DoseTargetErr = DoseTargetErr *doseConversion

layersize = 3
energy_notarget = 0
energy_heterotarget = 0
    
if(targetSelect == 0):
    print(f"Total energy deposition no target: {energy_notarget} MeV")
if(targetSelect == 1):
    print(f"Total energy deposition homogeneous target: {energy_notarget} MeV")
if(targetSelect == 2):
    print(f"Total energy deposition no target: {energy_notarget} MeV")
    print(f"Total energy deposition heterogeneous target: {energy_heterotarget} MeV")


plt.rcParams.update({'font.size': 20})
fig, ax1 = plt.subplots(figsize=(18, 13))
ax1.set_title(f'Fitted energy depth dose distribution {in_title[targetSelect]}')
if(bhetero):
    print(f'Fitted energy depth dose distribution {in_title[2]}')
    ax1.set_title(f'Fitted energy depth dose distribution with the heterogeneous target')
    ax1.set_title(f'Fitted energy depth dose distribution')
else:
    print(f'Fitted energy depth dose distribution {in_title[targetSelect]}')
    
start_time = time.time()

z = np.linspace(0, 40, 4001)

depth = np.array(depth[::-1])
Dose = Dose[::-1]
depthTarget=np.array(depthTarget[::-1])
DoseTarget= DoseTarget[::-1]
print(f"depth: {depth}")
print(f"Dose {Dose}")
print(f"depthTarget {depthTarget}")
print(f"DoseTarget {DoseTarget}")
# f = interp1d(depth, Dose, kind='linear', fill_value="extrapolate")
f = interp1d(depth, Dose, kind='cubic', fill_value="extrapolate")

popt, pcov =  curve_fit(lambda x, amp, mean, stddev: right_sided_convolution(f, lambda x2: gaussian(x2, amp, mean, stddev), x), depthTarget, DoseTarget, p0 = [1, 6, 0.3], bounds=((0.999, 0, 0), (1.0001, 10, 0.5)))

stddev = np.sqrt(np.diag(pcov))
    
t = popt[1]
o_t = stddev[1]   
    
sigmat = popt[2]
o_sigmat = stddev[2]

pmod = sigmat**2/t*10000
sigma_pmod = np.sqrt((2*sigmat/t**2*o_sigmat)**2+(sigmat**2/t**2*o_t)**2)*10000

print(f"amp = {popt[0]} +- {stddev[0]}")
print(f"t = {popt[1]} +- {stddev[1]}")
print(f"sigmat = {popt[2]} +- {stddev[2]}")
print(f"pmod = {pmod} +- {sigma_pmod}")

# popt = [1, 7.5, 0.3]
plt.errorbar(depth, Dose, DoseErr, depthErr, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color='#004600', label="No target data points") 
plt.errorbar(depthTarget, DoseTarget, DoseTargetErr, depthTargetErr, fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color="#cc7000", label="Hetero. data points")
plt.plot(z, f(z), label='Linear interpolation of no target data')
plt.plot(z, right_sided_convolution(f, lambda x2: gaussian(x2, *popt), z), label='Right sided convolution: \n' fr"$t= {t:.3f}" "\pm" fr"{o_t:.3f}~cm, \sigma={sigmat:.3f} \pm {o_sigmat:.3f}~cm,~P_{{mod}}={pmod:.3e} \pm {sigma_pmod:.3f}$" rf"$~\mu m$")
plt.grid(True)
plt.xlabel('Water Equivalent Depth / cm')
plt.ylabel('Energy Deposition per thickness  / MeV/mm')
plt.legend(loc='upper left',  fancybox=False, edgecolor='black')
fig.tight_layout()
plt.show()