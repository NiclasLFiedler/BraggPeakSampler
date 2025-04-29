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

def depth_dose_convolved(x, Phi0, R0, sigma, epsilon, amp, mean, stddev):
    return right_sided_convolution(lambda x1: depth_dose_distribution(x1, Phi0, R0, sigma, epsilon), lambda x2: gaussian(x2, amp, mean, stddev), x)


def convolution_fit(Xconv, Yconv, params, weights):
    popt, pcov = curve_fit(lambda x, amp, mean, stddev: depth_dose_convolved(x, params.Phi0, params.R0, params.sigma, params.epsilon, amp, mean, stddev), xdata= Xconv, ydata=Yconv, sigma=weights, p0 = [1, 7, 0.5], bounds=((0.9, 0, 0), (1.1, 10, 1)))
    errors = np.sqrt(np.diag(pcov))
    return [popt, errors]

def cyl_gauss(a,x):
    xarr=np.copy(x)
    yarr = np.copy(xarr)
    branch1 = -20.0
    branch2 = 20.0

    indset1 = (xarr<branch1) 
    x1 = xarr[indset1]
    y1 = np.sqrt(2*np.pi)/special.gamma(-a)*(-x1)**(-a-1)
    yarr[indset1] = y1

    indset2 = ((xarr>=branch1) & (xarr<branch2))
    x2 = xarr[indset2]
    y2a = special.pbdv(a,x2)[0]
    y2b = np.exp(-x2*x2/4)
    y2 = y2a*y2b
    yarr[indset2] = y2

    indset3 = (xarr>=branch2)
    yarr[indset3] = 0.0

    return yarr

def depth_dose_distribution(z, Phi0, R0, sigma, epsilon): 
    beta =  0.0
    gamma = 0.6

    return Phi0*sigma**(1/p_h2o)*special.gamma(1/p_h2o)/(np.sqrt(2*np.pi)*p_h2o*a_h2o**(1/p_h2o)*(1+beta*R0))*(sigma**(-1)*cyl_gauss(-1/p_h2o,(z-R0)/sigma)+(beta/p_h2o + gamma*beta + epsilon/R0)*cyl_gauss(-1/p_h2o-1,(z-R0)/sigma))

def depth_dose_distribution_ionization(z, Phi0, R0, sigma, epsilon): 
    beta =  0.0

    return Phi0*sigma**(1/p_h2o)*special.gamma(1/p_h2o)/(np.sqrt(2*np.pi)*p_h2o*a_h2o**(1/p_h2o)*(1+beta*R0))*(sigma**(-1)*cyl_gauss(-1/p_h2o,(z-R0)/sigma))

def depth_dose_distribution_nonelastic(z, Phi0, R0, sigma, epsilon): 
    beta =  0.0
    gamma = 0.6

    return Phi0*sigma**(1/p_h2o)*special.gamma(1/p_h2o)/(np.sqrt(2*np.pi)*p_h2o*a_h2o**(1/p_h2o)*(1+beta*R0))*((beta/p_h2o + gamma*beta + epsilon/R0)*cyl_gauss(-1/p_h2o-1,(z-R0)/sigma))

def gaussian_with_cutoff(mean, sigma, cutoff=2.5):
    while True:
        value = np.random.normal(mean, sigma)
        if mean-cutoff <= value <= mean+cutoff:
            return value

def range_energy_relationship(E, alpha, p):
    return alpha*E**p

def bortfeld_fit(x, y, Phi0, R0, epsilon, sigma, weights=None):    
    if weights is None:
        weights = np.ones_like(y)
    sigma_weights = 1 / weights
    
    params = fit_params()
    popt, pcov = curve_fit(lambda z, Phi0, R0, sigma, epsilon: depth_dose_distribution(z, Phi0, R0, sigma, epsilon), x, y, p0=[Phi0, R0, sigma*10, epsilon], bounds=((Phi0*0.5, R0 - 3*sigma, 0.1*sigma, 0), (Phi0*1.5, R0 + 3.5*sigma, 10*sigma, 1)), sigma=sigma_weights, maxfev=int(1e8))
    params.curve = depth_dose_distribution(z, *popt)
    params.Phi0 = popt[0] 
    params.R0 = popt[1] 
    params.sigma = popt[2]
    params.epsilon = popt[3]
    params.stddev = np.sqrt(np.diag(pcov))
    return params

def bortfeld_fit_hetero(x, y, Phi0, R0, sigma, epsilon, weights=None):    
    if weights is None:
        weights = np.ones_like(y)
    sigma_weights = 1 / weights
    
    params = fit_params()
    #popt, _ = curve_fit(lambda z, R0, sigma: depth_dose_distribution(z, Phi0, R0, sigma, epsilon), x, y, p0=[R0, sigma], bounds=((R0 - 3*sigma, 0.1*sigma), (R0 + 3.5*sigma, 3*sigma)), sigma=sigma_weights, maxfev=int(1e8))
    popt, pcov = curve_fit(lambda z, Phi0, R0, sigma: depth_dose_distribution(z, Phi0, R0, sigma, epsilon), x, y, p0=[Phi0, R0, sigma], bounds=((Phi0*0.5, R0 - 3*sigma, 0.01*sigma), (Phi0*1.5, R0 + 3.5*sigma, 3*sigma)), sigma=sigma_weights, maxfev=int(1e8))
    params.curve = depth_dose_distribution(z, popt[0], popt[1], popt[2], epsilon)
    params.Phi0 = popt[0] 
    params.R0 = popt[1] 
    params.sigma = popt[2]
    params.epsilon = epsilon
    params.stddev = np.sqrt(np.diag(pcov))
    return params

def characterize_z_D_curve(z, D):
    results = {}
    D100_index = np.argmax(D)
    D100       = D[D100_index]
    R100       = z[D100_index]
    results.update({'D100': D100, 'R100': R100})

    z_proximal    = z[:D100_index]
    dose_proximal = D[:D100_index]
    z_distal      = z[D100_index:]
    dose_distal   = D[D100_index:]

    R90P = z_proximal[np.argmin(np.abs(dose_proximal - 0.9 * D100))]
    R90D = z_distal  [np.argmin(np.abs(dose_distal   - 0.9 * D100))]
    R80P = z_proximal[np.argmin(np.abs(dose_proximal - 0.8 * D100))]
    R80D = z_distal  [np.argmin(np.abs(dose_distal   - 0.8 * D100))]
    R50P = z_proximal[np.argmin(np.abs(dose_proximal - 0.5 * D100))]
    R50D = z_distal  [np.argmin(np.abs(dose_distal   - 0.5 * D100))]
    R20D = z_distal  [np.argmin(np.abs(dose_distal   - 0.2 * D100))]
    R10D = z_distal  [np.argmin(np.abs(dose_distal   - 0.1 * D100))]
    results.update({'R90P': R90P, 'R90D': R90D, 'R80P': R80P, 'R80D': R80D, 'R50P': R50P, 'R50D': R50D, 'R20D': R20D, 'R10D': R10D})

    FWHM    = R50D  - R50P
    DFO2080 = R20D  - R80D
    DFO1090 = R10D  - R90D
    results.update({'FWHM': FWHM, 'DFO2080': DFO2080, 'DFO1090': DFO1090})

    return results

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
epsilon = 0.001*beamEnergy
Phi0 = max(y_data1)/20

beta = 0.012

resolution = 0.01*np.min(np.diff(z))

spline_func = interp1d(x_data, y_data1, kind='cubic')
z_spline    = np.linspace(min(x_data), max(x_data), round((max(x_data)-min(x_data)) / resolution ))
quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))


R0 = quantities['R80D']
if file == "notarget":
    R0 = range_energy_relationship(beamEnergy, a_h2o, p_h2o)

sigmaMono = (beta*R0**0.935)
sigmaE0   = 0.01*beamEnergy
sigma     = np.sqrt(sigmaMono**2+sigmaE0**2*a_h2o**2*p_h2o**2*beamEnergy**(2*p_h2o-2))

top_indices = np.argsort(y_data1)[-3:]
weights = np.ones_like(y_data1)
weights[top_indices] = 1


ax1.errorbar(x_data, y_data1, y_sigma1, x_sigma, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color='#004600', label="No target data points") 
ax1.errorbar(x_data2, y_data2, y_sigma2, x_sigma2, fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color="#cc7000", label="Hetero. data points") #Convolution

# Convolve
x = np.asarray(x_data)
y2 = gaussian(x, 1, 2.55, 0.277)
print(y2)
smoothed = right_sided_convolution(interp1d(x_data, y_data1, kind='linear', fill_value="extrapolate"), interp1d(x_data, y2, kind='linear', fill_value="extrapolate"), x)

# Plot
plt.plot(x, smoothed, label='Right-sided convolved (Gaussian)')
plt.legend()
plt.show()