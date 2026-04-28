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
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import Akima1DInterpolator
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

def normalized_gaussian(x, amp, mean, stddev):
    """Properly normalized — integrates to 1. No amp parameter."""
    return amp * np.exp(-(x - mean)**2 / (2 * stddev**2)) / (np.sqrt(2 * np.pi) * stddev)

def right_sided_convolution_stable(f, amp, t_shift, sigma, z_values):
    """
    Convolves f(z) with a normalized Gaussian g(t) = N(t_shift, sigma),
    integrated over t in [0, t_max] where t_max adapts to sigma.
    
    Physical meaning: f is the no-target depth-dose, g describes
    the straggling + mean shift from the heterogeneous target.
    """
    z_values = np.asarray(z_values)
    
    # Adaptive integration range: cover at least 5 sigma beyond mean
    t_max = t_shift + 6 * sigma
    t_max = max(t_max, 2.0)  # minimum range
    n_points = max(1000, int(t_max / sigma * 200))  # density scales with sigma
    
    t = np.linspace(0, t_max, n_points)
    g_vals = normalized_gaussian(t, amp,t_shift, sigma)
    
    # Vectorized: shape (len(z_values), len(t))
    Z, T = np.meshgrid(z_values, t, indexing='ij')
    integrand = f(Z + T) * g_vals  # g_vals broadcasts over z axis
    
    result = simpson(integrand, x=t, axis=1)
    return result

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

def fft_convolution_onesided(f_interp, amp, t_shift, sigma, z_values):
    z_values = np.asarray(z_values)
    
    z_min = z_values.min()
    z_max = z_values.max()
    dz = 0.01
    # dz = z_values[1]-z_values[0]
    z_grid = np.arange(z_min, z_max, dz)
    
    f_vals = f_interp(z_grid)
    # plt.plot(z_grid, f_vals, label='f(z) on grid')
    # plt.show()
    # Kernel on positive t axis — this shifts f to the LEFT when convolved
    t_kernel = np.arange(0, t_shift + 8 * sigma, dz)
    kernel = np.exp(-(t_kernel - t_shift)**2 / (2 * sigma**2))
    kernel /= kernel.sum() * dz
    kernel *= amp
    # Use 'valid'-style: we want (f * g)[z] = sum_t f(z+t)*g(t)
    # fftconvolve with mode='full' gives sum_t f(z-t)*g(t), i.e. wrong direction
    # Flip f to get the correct left-shift direction
    conv = fftconvolve(f_vals[::-1], kernel * dz, mode='full')
    conv = conv[::-1]  # flip back
    
    z_conv = np.arange(len(conv)) * dz + z_grid[0] - (len(kernel) - 1) * dz
    
    conv_interp = interp1d(z_conv, conv, kind='linear',
                           bounds_error=False, fill_value=0.0)
    return conv_interp(z_values)

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
targetColorMap = ["#1f77b4", "#4e79a7", "#76b7b2", "#bab0ac", "#f28e2b", "#e15759", "#9c755f"]

dataFile = np.load(f"{dataset}/{file}/input/depthdose{nLayers}.npz")
print(f"File: {dataFile}")

depth = dataFile["depth"]
depthErr = dataFile["depth_err"]
Dose = dataFile["dose"]
DoseErr = dataFile["dose_err"]
unfoldedEntries = dataFile["dose"][0]

targetThicknesses = [200]
resolution = []

if(bhetero):
    for thickness in targetThicknesses:
        dataFile = np.load(f"{dataset}/{targetFile}/input/depthdose{nLayers}_{thickness}.npz")
        print(f"Hetero File: {dataFile}")
        depthTarget.append(dataFile["depth"])
        depthTargetErr.append(dataFile["depth_err"])
        DoseTarget.append(dataFile["dose"])
        DoseTargetErr.append(dataFile["dose_err"])
        unfoldedEntriesTarget.append(dataFile["dose"][0])

print(depthTarget)
print(depthTargetErr)
print(DoseTarget)
print(DoseTargetErr)

if(targetSelect == 1):
    for thickness in targetThicknesses:
        dataFile = np.load(f"{dataset}/{targetFile}/input/depthdose{nLayers}_{targetThickness}.npz")
        depthTarget.append(dataFile["depth"])
        depthTargetErr.append(dataFile["depth_err"])
        DoseTarget.append(dataFile["dose"])
        DoseTargetErr.append(dataFile["dose_err"])
        unfoldedEntriesTarget.append(dataFile["dose"][0])


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

z = np.linspace(0, crystalSize[2], 4001)

f_pchip = interp1d(depth, Dose, kind='linear', fill_value="extrapolate")
f_pchip = interp1d(depth, Dose, kind='cubic', fill_value="extrapolate")
# f_pchip = PchipInterpolator(depth, Dose, extrapolate=False)
# f_pchip = Akima1DInterpolator(depth, Dose)

def interpolate(z):
    z = np.asarray(z)
    result = f_pchip(z)
    result = np.where(np.isnan(result), 0.0, result)  # zero outside data range
    result = np.clip(result, 0.0, None)               # no negative values
    return result



# f = interpolate(z)
# # popt, pcov =  curve_fit(lambda x, amp, mean, stddev: right_sided_convolution(f, lambda x2: gaussian(x2, amp, mean, stddev), x), depthTarget, DoseTarget, p0 = [1, 6, 0.3], bounds=((0.999, 0, 0), (10, 10, 0.5)), maxfev=100000)

fitParamsRaw = []
fitParams = []

for index, thickness in enumerate(targetThicknesses):
    fitParams.append(curve_fit(
        lambda x, amp, t, sigma: fft_convolution_onesided(interpolate, amp, t, sigma, x),
        depthTarget[index],
        DoseTarget[index],
        p0=[1, 6, 0.03],
        bounds=((0.999, 0, 1e-10), (2, 10, 1)),
        maxfev=100000,
    ))

    stddev = np.sqrt(np.diag(fitParams[index][1]))
       
    t = fitParams[index][0][1]
    o_t = stddev[1]   

    sigmat = fitParams[index][0][2]
    o_sigmat = stddev[2]

    pmod = sigmat**2/t*10000
    sigma_pmod = np.sqrt((2*sigmat/t**2*o_sigmat)**2+(sigmat**2/t**2*o_t)**2)*10000

    # print(f"Target Thickness {thickness}")
    # print(f"amp = {fitParams[index][0][0]} +- {stddev[0]}")
    # print(f"t = {t} +- {stddev[1]}")
    print(f"Thickness {thickness}: sigmat = {sigmat} +- {stddev[2]}")
    # print(f"pmod = {pmod} +- {sigma_pmod}")

# popt = [1, 7.5, 0.3]
plt.errorbar(depth, Dose, DoseErr, depthErr, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color=targetColorMap[0], label="No target data points") 
plt.plot(z, interpolate(z), label='Linear interpolation of no target data')

for index, thickness in enumerate(targetThicknesses):
    plt.errorbar(depthTarget[index], DoseTarget[index], DoseTargetErr[index], depthTargetErr[index], fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color=targetColorMap[index+1], label="Hetero. data points")
    plt.plot(z, fft_convolution_onesided(interpolate, *fitParams[index][0], z),color=targetColorMap[index+1])#, label='Right sided convolution: \n' fr"$t= {t:.3f}" r"\pm" fr"{o_t:.3f}~cm, \sigma={sigmat:.3f} " r"\pm" rf" {o_sigmat:.3f}~cm,~P_{{mod}}={pmod:.3e}" r"\pm" fr"{sigma_pmod:.3f}$" rf"$~\mu m$")

plt.grid(True)
plt.xlabel('Water Equivalent Depth / cm')
plt.ylabel('Energy Deposition per thickness  / MeV/mm')
plt.legend(loc='upper left',  fancybox=False, edgecolor='black')
fig.tight_layout()
plt.show()