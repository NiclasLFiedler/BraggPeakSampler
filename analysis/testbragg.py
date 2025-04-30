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
            zlin = np.linspace(0,40,400)
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

with open("config.json", "r") as file:
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
discardIndex = config["discardIndex"]

datasets = ["MIT_05_2024", "simulation"]
in_data = ["notarget", "homotarget", "heterotarget"]
in_title = ["without a target", "with the homogeneous target", "with the heterogeneous target"]

dataset = datasets[datasetSelect]
file = in_data[targetSelect]
if targetSelect == 2: file = in_data[0]

nosave = False

bhetero = False
if targetSelect == 2: bhetero = True
nbOfFits = 0
nbOfFitsHetero = 0
lineWidth = 2
capSize = 3

targetfile = uproot.open(f"../data/{dataset}/{file}/output/{file}Means.root")
targettree = targetfile["meantree"]
y_data1 = targettree["mean"].array().to_numpy()
y_sigma1 = targettree["error"].array().to_numpy()
x_data = targettree["x"].array().to_numpy()
x_data = [x/10 for x in x_data]
x_sigma = targettree["x_sigma"].array().to_numpy()
x_sigma = [x/10 for x in x_sigma]

if(bhetero):
    hetero_file = in_data[2]
    heterotargetfile = uproot.open(f"../data/{dataset}/{hetero_file}/output/{hetero_file}Means.root")
    heterotree = heterotargetfile["meantree"]
    y_data2 = heterotree["mean"].array().to_numpy()
    y_sigma2 = heterotree["error"].array().to_numpy()
    x_data2 = heterotree["x"].array().to_numpy()
    x_data2 = [x/10 for x in x_data2]
    x_sigma2 = heterotree["x_sigma"].array().to_numpy()
    x_sigma2 = [x/10 for x in x_sigma2]

#plt.style.use(['science','notebook','grid'])
plt.rcParams.update({'font.size': 20})
fig, ax1 = plt.subplots(figsize=(16, 8))
ax1.set_title(f'Fitted energy depth dose distribution {in_title[targetSelect]}')
if(bhetero):
    print(f'Fitted energy depth dose distribution {in_title[2]}')
    ax1.set_title(f'Fitted energy depth dose distribution with the heterogeneous target')
    ax1.set_title(f'Fitted energy depth dose distribution')
else:
    print(f'Fitted energy depth dose distribution {in_title[targetSelect]}')
    
start_time = time.time()

z = np.linspace(0, 40, 4001)

ax1.errorbar(x_data, y_data1, y_sigma1, x_sigma, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color='#004600', label="No target data points") 
ax1.errorbar(x_data2, y_data2, y_sigma2, x_sigma2, fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color="#cc7000", label="Hetero. data points") #Convolution

f = interp1d(x_data, y_data1, kind='linear', fill_value="extrapolate")

popt, pcov =  curve_fit(lambda x, amp, mean, stddev: right_sided_convolution(f, lambda x2: gaussian(x2, amp, mean, stddev), x), x_data2, y_data2, p0 = [1, 2, 0.2], bounds=((0.9, 0, 0), (1.1, 10, 1)))

print(f"amp = {popt[0]}")
print(f"t = {popt[1]}")
print(f"sigmat = {popt[2]}")
print(f"pmod = {popt[2]**2/popt[1]*10**4}")
plt.plot(z, f(z), label='Right-sided convolved (Gaussian)')
plt.plot(z, right_sided_convolution(f, lambda x2: gaussian(x2, *popt), z), label='Right-sided convolved (Gaussian)')
plt.legend()
plt.show()