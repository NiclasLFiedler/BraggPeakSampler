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

datasets = ["MIT_05_2024", "simulation", "beamtime"]
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

y_data1 = [5.76817, 5.78830, 5.91962, 6.02027, 6.17647, 6.31865, 6.43451, 6.57329, 6.79197, 7.03662, 7.28912, 4.70684, 4.97107, 5.04074, 5.33586, 5.69500, 5.73153, 5.96638, 6.15702, 6.33759, 6.89649, 7.27459, 7.80804, 8.40822, 9.70558, 10.8170, 13.4784, 20.2835, 11.1428, 5.36193, 5.45272, 5.38791]
y_sigma1 = [1.05952, 1.00481, 1.04628, 1.05401, 1.17253, 1.02417, 1.19482, 1.08550, 1.34963, 1.22602, 1.41804, 1.75369, 1.83540, 2.31214, 1.88457, 1.50985, 2.35451, 3.49368, 1.52134, 2.08780, 1.82112, 1.73598, 2.25831, 2.29804,  4.19429, 2.48435, 2.27875, 6.30751, 2.37186, 1.24200, 1.17689, 9.92402]
counts = [197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321,197321, 197321, 197321, 197321, 10486, 1546, 1541, 763] 

if(targetSelect == 1):
    y_data1 = [6.00176, 6.07318, 6.39285, 6.43914, 6.71883, 6.45248, 7.03723, 7.33222, 7.58636, 8.02201, 8.25445, 6.20799, 5.40542, 5.38858, 6.25463, 6.75024, 5.80082, 6.20833, 9.20776, 10.8855, 15.5453, 9.7055, 2.05837, 7.7235, 8.14947, 13.885, 13.5339, 18.5561, 12.1416, 4.57177, 4.13783, 4.73284]
    y_sigma1 = [1.08812, 1.05688, 1.05828, 1.1409, 1.15742, 0.931626, 1.24063, 0.943488, 1.42111, 1.27323, 1.50352, 2.6597, 1.83007, 2.23962, 1.23171, 1.43868, 2.28398, 3.55953, 1.93801, 2.78693, 4.27617, 2.19289, 1.43113, 2.76886, 4.17295, 5.40727, 5.76706, 6.11727, 5.10709, 1.14239, 0.986111, 0.949723]
    counts = [300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 300000, 95108, 5253, 2939, 16824, 328, 125, 387, 664, 616, 730, 294] 

y_data2 = [5.98458, 6.03263, 6.36008, 6.38608, 6.68541, 6.44159, 7.00916, 7.27892, 7.54899, 7.98059, 8.0964, 4.49038, 5.3649, 5.32143, 6.21623, 6.69112, 5.76967, 6.75673, 8.76714, 9.60902, 13.4133,12.8486, 2.03278, 8.096, 12.5796, 13.8951, 13.2846, 18.0847, 10.8717, 4.69535, 4.3318, 6.49868]
y_sigma2 = [1.09676, 1.04085, 1.07188, 1.08782, 1.11803, 0.859022, 1.24689, 0.97747, 1.39625, 1.25349, 1.50665, 2.5856, 1.8828, 2.23713, 1.1592, 1.45099, 2.4415, 3.7286, 1.91762, 2.90847, 4.83249, 5.40193, 1.44396, 2.90994, 13.4627, 4.41395, 2.78418, 7.49992, 3.67727, 1.21208, 1.03319, 1.39521]
counts2 = [250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 250538, 167869, 31879, 2937, 295, 217, 143, 369, 800, 736, 897, 383] 

layersize = 3
energy_notarget = 0
energy_heterotarget = 0

for index, value in enumerate(y_data1):
    if(index == 11):
        layersize = 2
    y_data1[index] = value*counts[index]/counts[0]*1/layersize
    energy_notarget += y_data1[index]*layersize
    y_sigma1[index] = y_sigma1[index]*counts[index]/counts[0]*1/layersize
    if(targetSelect == 2):  
        y_data2[index] = y_data2[index]*counts2[index]/counts2[0]*1/layersize
        energy_heterotarget += y_data2[index]*layersize
        y_sigma2[index] = y_sigma2[index]*counts2[index]/counts2[0]*1/layersize
    
if(targetSelect == 0):
    print(f"Total energy deposition no target: {energy_notarget} MeV")
if(targetSelect == 1):
    print(f"Total energy deposition homogeneous target: {energy_notarget} MeV")
if(targetSelect == 2):
    print(f"Total energy deposition no target: {energy_notarget} MeV")
    print(f"Total energy deposition heterogeneous target: {energy_heterotarget} MeV")

# for _ in range(20):
#     y_data1 = np.append(y_data1, 0)
#     y_sigma1 = np.append(y_sigma1, 0)
#     new_x = x_data[-1] + x_data[-1] - x_data[-2]
#     x_data = np.append(x_data, new_x)
#     x_sigma = np.append(x_sigma, 0)

#     y_data2 = np.append(y_data2, 0)
#     y_sigma2 = np.append(y_sigma2, 0)
#     new_x = x_data2[-1] + x_data2[-1] - x_data2[-2]
#     x_data2 = np.append(x_data2, new_x)
#     x_sigma2 = np.append(x_sigma2, 0)

#plt.style.use(['science','notebook','grid'])
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

ax1.errorbar(x_data, y_data1, y_sigma1, x_sigma, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color='#004600', label="No target data points") 
ax1.errorbar(x_data2, y_data2, y_sigma2, x_sigma2, fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color="#cc7000", label="Hetero. data points") #Convolution

f = interp1d(x_data, y_data1, kind='linear', fill_value="extrapolate")
# f = interp1d(x_data, y_data1, kind='cubic', fill_value="extrapolate")

popt, pcov =  curve_fit(lambda x, amp, mean, stddev: right_sided_convolution(f, lambda x2: gaussian(x2, amp, mean, stddev), x), x_data2, y_data2, p0 = [1, 2, 0.2], bounds=((0.8, 0, 0), (1.1, 10, 1)))

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

plt.plot(z, f(z), label='Linear interpolation of no target data')
plt.plot(z, right_sided_convolution(f, lambda x2: gaussian(x2, *popt), z), label='Right sided convolution: \n' fr"$t= {t:.3f}" "\pm" fr"{o_t:.3f}~cm, \sigma={sigmat:.3f} \pm {o_sigmat:.3f}~cm,~P_{{mod}}={pmod:.3f} \pm {sigma_pmod:.3f}$" rf"$~\mu m$")
plt.grid(True)
plt.xlabel('Water Equivalent Depth / cm')
plt.ylabel('Energy Deposition per thickness  / MeV/mm')
plt.legend(loc='upper left',  fancybox=False, edgecolor='black')
fig.tight_layout()
plt.show()