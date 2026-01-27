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
    popt, pcov = curve_fit(lambda x, amp, mean, stddev: depth_dose_convolved(x, params.Phi0, params.R0, params.sigma, params.epsilon, amp, mean, stddev), xdata= Xconv, ydata=Yconv, sigma=weights, p0 = [1, 2, 0.2], bounds=((0.9, 0, 0), (1.1, 10, 2)))
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

def bortfeld_fit(x, y, Phi0, R0, Sigma, epsilon, weights=None):    
    if weights is None:
        weights = np.ones_like(y)
    sigma_weights = 1 / weights
    
    params = fit_params()
    popt, pcov = curve_fit(lambda z, Phi0, R0, Sigma, epsilon: depth_dose_distribution(z, Phi0, R0, Sigma, epsilon), x, y, p0=[Phi0, R0, Sigma*10, 0.1*epsilon], bounds=((Phi0*0.25, R0 - 3*Sigma, 0.001*Sigma, 0), (Phi0*1.5, R0 + 3.5*Sigma, 10*Sigma, 2)), sigma=sigma_weights, maxfev=int(1e8))
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
targetThickness = config["targetThickness"]

datasets = ["MIT_05_2024", "simulation", "paperBeamtime"]
in_data = ["notarget", "homotarget", "heterotarget"]
in_title = ["without a target", "with the homogeneous target", "with the heterogeneous target"]

dataset = datasets[datasetSelect]
file = in_data[0]
targetFile = in_data[targetSelect]
# if(targetSelect == 1):
#     targetFile = f"{targetFile}{targetThickness}"
nosave = False

bhetero = False
if targetSelect == 2: bhetero = True

nbOfFits = 0
nbOfFitsHetero = 0
lineWidth = 2
capSize = 3

depth = []
depthTarget = []
depthTargetTemp = []
depthErr = []
depthErrTarget = []
depthErrTargetTemp = []
dose = []
doseTarget = []
doseTargetTemp = []
doseErr = []
doseErrTarget = []
doseErrTargetTemp = []
unfoldedDose = []
unfoldedDoseTarget = []
unfoldedDoseTargetTemp = []
unfoldedDoseErr = []
unfoldedDoseErrTarget = []
unfoldedDoseErrTargetTemp = []

unfoldedEntries = []
unfoldedEntriesTarget = []
unfoldedEntriesTargetTemp = []

WETConv = 4.73399
WETConvTeflon = 0.289088
WETConvAlu = 0.289088
t_n = 0
x = 3
xTeflon = 0.5
xAlu = 0.02
WETDepth = []
WETDepthErr = []
WETDepthErrTarget = []
# for ch in range(nLayers):
#     if ch == 11:
#         x = 2
#     WETDepth.append(t_n+x/20*WETConv+xTeflon/20*WETConvTeflon+xAlu/20*WETConvAlu)
#     t_n += x*WETConv/10+xTeflon/20*WETConvTeflon+xAlu/20*WETConvAlu0
#     WETDepthErr.append(x/10*WETConv/np.sqrt(12))

WETDepth = [7.155416420506668, 21.455056176722735, 35.73121953386203, 50.00837171548292, 64.28558352543337, 78.54541699109387, 92.7836041434076, 107.00308108188163, 121.21292882757369, 135.4157474876818, 149.61160933141255, 161.48639261416054, 171.05138584768923, 180.61204957657964, 190.15178624767393, 199.66085509299435, 209.1565117771304, 218.64418551465582, 228.12356365600664, 237.59409600831475, 247.03996135329004, 256.4575830284209, 265.84166488741903, 275.19102039365407, 284.4969697433472, 293.7399813658173, 302.7956853321799, 311.9429171602238, 321.67214994835393, 331.6915553001954, 341.61841079794243, 351.4893609799391]

WET = [14.310832841013337, 28.599279512432133, 42.86315955529193, 57.15358387567392, 71.41758317519282, 85.67325080699493, 99.89395747982027, 114.11220468394299, 128.3136529712044, 142.51784200415926, 156.70537665866584, 166.26740856965526, 175.83536312572323, 185.38873602743604, 194.9148364679118, 204.4068737180769, 213.90614983618389, 223.38222119312775, 232.86490611888553, 242.32328589774397, 251.7566368088361, 261.1585292480057, 270.5248005268324, 279.8572402604758, 289.1366992262187, 298.3432635054159, 307.24810715894387, 316.6377271615037, 326.7065727352042, 336.67653786518656, 346.5602837306983, 356.41843822917986]

WETTarget = [np.float64(14.260531500089058), np.float64(28.46471731352301), np.float64(42.664750730649416), np.float64(56.832721552397004), np.float64(70.98226539002906), np.float64(85.16940527757126), np.float64(99.34093980739712), np.float64(113.53812826046382), np.float64(127.66679859692594), np.float64(141.80430163730276), np.float64(155.90094177435978), np.float64(165.3564491136862), np.float64(174.8355735832181), np.float64(184.268456238946), np.float64(193.67271340024672), np.float64(203.06511871147984), np.float64(212.39759846374548), np.float64(221.7005628885991), np.float64(230.94648960872445), np.float64(240.1028519371951), np.float64(249.05982850300109), np.float64(258.59018765964584), np.float64(268.5356069368744), np.float64(278.4626318829971), np.float64(288.3119606056622), np.float64(298.17920282074414), np.float64(308.0464450358261), np.float64(317.913687250908), np.float64(327.78092946598997), np.float64(337.6481716810719), np.float64(347.51541389615386), np.float64(357.3826561112358)]

WETDepthTarget = [np.float64(7.130265750044529), np.float64(21.362624406806034), np.float64(35.56473402208621), np.float64(49.74873614152321), np.float64(63.90749347121303), np.float64(78.07583533380016), np.float64(92.25517254248419), np.float64(106.43953403393047), np.float64(120.60246342869488), np.float64(134.73555011711434), np.float64(148.85262170583127), np.float64(160.62869544402298), np.float64(170.09601134845215), np.float64(179.55201491108204), np.float64(188.97058481959635), np.float64(198.36891605586328), np.float64(207.73135858761265), np.float64(217.0490806761723), np.float64(226.32352624866178), np.float64(235.52467077295978), np.float64(244.5813402200981), np.float64(253.82500808132346), np.float64(263.5628972982601), np.float64(273.4991194099357), np.float64(283.3872962443296), np.float64(293.24558171320314), np.float64(303.11282392828514), np.float64(312.980066143367), np.float64(322.847308358449), np.float64(332.7145505735309), np.float64(342.5817927886129), np.float64(352.4490350036948)]

for ch in range(nLayers):
    WETDepth[ch] = WETDepth[ch]/10 
    if ch == 0:
        WETDepthErr.append((WET[ch])/(np.sqrt(12)*10))
    else:
        WETDepthErr.append((WET[ch]-WET[ch-1])/(np.sqrt(12)*10))

for ch in range(nLayers):
    WETDepthTarget[ch] = WETDepthTarget[ch]/10 
    if ch == 0:
        WETDepthErrTarget.append((WETTarget[ch])/(np.sqrt(12)*10))
    else:
        WETDepthErrTarget.append((WETTarget[ch]-WETTarget[ch-1])/(np.sqrt(12)*10))

for ch in range(nLayers):
    dataFile = np.load(f"../data/{dataset}/{file}/output/unfold/data/Unfolded{ch}.npz")
    # Tdepth, Tdeptherr = dataFile["depth"] 
    # depth.append(Tdepth)
    # depthErr.append(Tdeptherr)
    depth.append(WETDepth[ch])
    depthErr.append(WETDepthErr[ch])
    
    x_pred, y_pred = dataFile["predicted_measured"]
    x_meas, y_meas = dataFile["measured_energy"]
    mask = (y_meas > 0) & (y_pred > 0)

    ndof = len(y_meas[mask])  # no fitted params
    sigma = np.sqrt(y_meas[mask])

    chi2 = np.sum((y_pred[mask] - y_meas[mask])**2 / sigma**2)
    chi2_red = chi2 / ndof

    print(f"Ch: {ch} Reduced Chi^2: {chi2_red}")
    Tdose, Tdoseerr = dataFile["unfoldedDose"]
    unfoldedDose.append(Tdose)
    unfoldedDoseErr.append(Tdoseerr)
    Tdose, Tdoseerr = dataFile["dose"]
    dose.append(Tdose)
    doseErr.append(Tdoseerr)
    unfoldedEntries.append(np.sum(dataFile["unfolded"][1]))

if(bhetero):
    for ch in range(nLayers):
        dataFile = np.load(f"../data/{dataset}/{targetFile}/output/unfold/data/Unfolded{ch}.npz")

        # Tdepth, Tdeptherr = dataFile["depth"] 
        # depthTarget.append(Tdepth)
        # depthErrTarget.append(Tdeptherr)
        
        depthTarget.append(WETDepthTarget[ch])
        depthErrTarget.append(WETDepthErrTarget[ch])
        Tdose, Tdoseerr = dataFile["unfoldedDose"]
        unfoldedDoseTarget.append(Tdose)
        unfoldedDoseErrTarget.append(Tdoseerr)
        Tdose, Tdoseerr = dataFile["dose"]
        doseTarget.append(Tdose)
        doseErrTarget.append(Tdoseerr)
        unfoldedEntriesTarget.append(np.sum(dataFile["unfolded"][1]))

targetThicknesses = [52, 53, 54, 55, 56, 57]
if(targetSelect == 1):
    for thickness in targetThicknesses:
        depthTargetTemp = []
        depthErrTargetTemp = []
        unfoldedDoseErrTargetTemp = []
        unfoldedDoseTargetTemp = []
        doseTargetTemp = []
        doseErrTargetTemp = []
        unfoldedEntriesTargetTemp = []
        for ch in range(nLayers):
            dataFile = np.load(f"../data/{dataset}/{targetFile}{thickness}/output/unfold/data/Unfolded{ch}.npz")

            # Tdepth, Tdeptherr = dataFile["depth"] 
            # depthTargetTemp.append(Tdepth)
            # depthErrTargetTemp.append(Tdeptherr)
            
            depthTargetTemp.append(WETDepthTarget[ch])
            depthErrTargetTemp.append(WETDepthErrTarget[ch])
            
            Tdose, Tdoseerr = dataFile["unfoldedDose"]
            unfoldedDoseTargetTemp.append(Tdose)
            unfoldedDoseErrTargetTemp.append(Tdoseerr)
            Tdose, Tdoseerr = dataFile["dose"]
            doseTargetTemp.append(Tdose)
            doseErrTargetTemp.append(Tdoseerr)
            unfoldedEntriesTargetTemp.append(np.sum(dataFile["unfolded"][1]))
            # if thickness == 54:
                # chi2 = dataFile["reduced_chi2"]
                # print(f"{ch} Chi^2: {chi2}")
        depthTarget.append(depthTargetTemp)
        depthErrTarget.append(depthErrTargetTemp)
        unfoldedDoseTarget.append(unfoldedDoseTargetTemp)
        unfoldedDoseErrTarget.append(unfoldedDoseErrTargetTemp)
        doseTarget.append(doseTargetTemp)
        doseErrTarget.append(doseErrTargetTemp)
        unfoldedEntriesTarget.append(unfoldedEntriesTargetTemp)

# print(unfoldedEntriesTarget[1])
if targetSelect == 1:
    targetFile = f"{targetFile}{52}"

plt.rcParams.update({'font.size': 32})
fig, ax1 = plt.subplots(figsize=(16, 12))
# ax1.set_title(f'Fitted energy depth dose distribution {in_title[targetSelect]}')
# if(bhetero):
#     print(f'Fitted energy depth dose distribution {in_title[2]}')
#     ax1.set_title(f'LN300 Phantom: Fitted energy depth dose distribution')
#     # ax1.set_title(f'Fitted energy depth dose distribution')
# elif(targetSelect == 1):
#     ax1.set_title(f'Fitted energy depth dose distribution with PMMA Phantoms')
# else:
#     print(f'Fitted energy depth dose distribution {in_title[targetSelect]}')
    
start_time = time.time()

z = np.linspace(0, 40, 4001)
epsilon = 0.001*beamEnergy
Phi0 = max(unfoldedDose)/20
beta = 0.012

resolution = 0.01*np.min(np.diff(z))

spline_func = interp1d(depth, unfoldedDose, kind='cubic')
z_spline    = np.linspace(min(depth), max(depth), round((max(depth)-min(depth)) / resolution ))
quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))


R0 = quantities['R80D']
#if file == "notarget":
#    R0 = range_energy_relationship(beamEnergy, a_h2o, p_h2o)
print("Expected range from R80D fit: ", range_energy_relationship(225, a_h2o, p_h2o))
print("Expected range from R80D fit: ", range_energy_relationship(235, a_pwo, p_pwo))
sigmaMono = (beta*R0**0.935)
sigmaE0   = 0.01*beamEnergy
sigma     = np.sqrt(sigmaMono**2+sigmaE0**2*a_h2o**2*p_h2o**2*beamEnergy**(2*p_h2o-2))

top_indices = np.argsort(unfoldedDose)[-3:]
weights = np.ones_like(unfoldedDose)
weights[top_indices] = 1

params_list = []
for i in range(nbOfFits):
    if i%100 == 0:
        print(f"Fit: {i}/{nbOfFits}")
    x_with_err = []
    y_with_err = []
    for index, mean in enumerate(depth):
        x_with_err.append(gaussian_with_cutoff(mean, depthErr[index]))
    for index, mean in enumerate(unfoldedDose):
        y_with_err.append(gaussian_with_cutoff(mean, unfoldedDoseErr[index]))
    params = bortfeld_fit(x_with_err, y_with_err, Phi0, R0, sigma, epsilon, weights)
    
    if(params.Phi0 == 0):
        continue
    else:
        params_list.append(params)
        params_list[-1].curve = depth_dose_distribution(z, params.Phi0, params.R0, params.sigma, params.epsilon)
        
bestfit_params = bortfeld_fit(depth, unfoldedDose, Phi0, R0, sigma, epsilon, weights)
if(targetSelect == 2):
    doseConversion = unfoldedEntries[0]/unfoldedEntriesTarget[0]*(1+beta*(bestfit_params.R0-200*0.029))/(1+beta*bestfit_params.R0)
    print(doseConversion)
    
    # unfoldedDoseTarget = np.array(unfoldedDoseTarget)
    # unfoldedDoseErrTarget = np.array(unfoldedDoseErrTarget)
    for idx, value in enumerate(unfoldedDoseTarget):
        unfoldedDoseTarget[idx] = unfoldedDoseTarget[idx] * doseConversion
        unfoldedDoseErrTarget[idx] = unfoldedDoseErrTarget[idx] *doseConversion

    # unfoldedDoseTarget = unfoldedDoseTarget * doseConversion
    # unfoldedDoseErrTarget = unfoldedDoseErrTarget *doseConversion
elif(targetSelect == 1):
    for idx, unfoldedEntriesTargetSingle in enumerate(unfoldedEntriesTarget):
        doseConversion = unfoldedEntries[0]/unfoldedEntriesTargetSingle[0]*(1+beta*(bestfit_params.R0-targetThickness*0.118))/(1+beta*bestfit_params.R0)
        print(doseConversion)
        
        # unfoldedDoseTargetTemp = np.array(unfoldedDoseTarget[idx])
        # unfoldedDoseErrTargetTemp = np.array(unfoldedDoseErrTarget[idx])
        for index, value in enumerate(unfoldedDoseTarget[idx]):
            unfoldedDoseTarget[idx][index] = unfoldedDoseTarget[idx][index] * doseConversion
            unfoldedDoseErrTarget[idx][index] = unfoldedDoseErrTarget[idx][index]*doseConversion
        # unfoldedDoseTarget[idx] = unfoldedDoseTargetTemp * doseConversion
        # unfoldedDoseErrTarget[idx] = unfoldedDoseErrTargetTemp *doseConversion

targetColorMap = ["#1f77b4", "#4e79a7", "#76b7b2", "#bab0ac", "#f28e2b", "#e15759", "#9c755f"]

if(bhetero):
    ax1.errorbar(depth, unfoldedDose, unfoldedDoseErr, depthErr, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color="#000000")#, label="No target data") 
    ax1.errorbar(depthTarget, unfoldedDoseTarget, unfoldedDoseErrTarget, depthErrTarget, fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color=targetColorMap[0])#, label="Hetero. data") 
    
    depthTarget =  depthTarget[:8] + depthTarget[9:]
    unfoldedDoseTarget =  unfoldedDoseTarget[:8] + unfoldedDoseTarget[9:]
    unfoldedDoseErrTarget =  unfoldedDoseErrTarget[:8] + unfoldedDoseErrTarget[9:]
    depthErrTarget = depthErrTarget[:8] + depthErrTarget[9:]
    
elif(targetSelect == 1):
    ax1.errorbar(depth, unfoldedDose, unfoldedDoseErr, depthErr, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color='#000000')#, label="No target data")
    for idx, thickness in enumerate(targetThicknesses):
        ax1.errorbar(depthTarget[idx], unfoldedDoseTarget[idx], unfoldedDoseErrTarget[idx], depthErrTarget[idx], fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color=targetColorMap[idx])#, label=f"{thickness} mm PMMA data")

        depthTarget[idx] =  depthTarget[idx][:8] + depthTarget[idx][9:]
        unfoldedDoseTarget[idx] =  unfoldedDoseTarget[idx][:8] + unfoldedDoseTarget[idx][9:]
        unfoldedDoseErrTarget[idx] =  unfoldedDoseErrTarget[idx][:8] + unfoldedDoseErrTarget[idx][9:]
        depthErrTarget[idx] = depthErrTarget[idx][:8] + depthErrTarget[idx][9:]
else:
    ax1.errorbar(depth, unfoldedDose, unfoldedDoseErr, depthErr, fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color="#000000")#, label="No target data") 

print(type(unfoldedDose))
print(type(unfoldedDoseTarget))

depth =  depth[:8] + depth[9:]
unfoldedDose =  unfoldedDose[:8] + unfoldedDose[9:]
unfoldedDoseErr =  unfoldedDoseErr[:8] + unfoldedDoseErr[9:]
depthErr = depthErr[:8] + depthErr[9:]

print(f"Phi0: {bestfit_params.Phi0} \pm {bestfit_params.stddev[0]}")
print(f"R0: {bestfit_params.R0} \pm {bestfit_params.stddev[1]}")
print(f"sigma: {bestfit_params.sigma} \pm {bestfit_params.stddev[2]}")
print(f"epsilon: {bestfit_params.epsilon} \pm {bestfit_params.stddev[3]}")

if(targetSelect > 0):
    ax1.plot(z, bestfit_params.curve, color="black", linewidth = lineWidth, label=fr"No target")# rf"$R_0={bestfit_params.R0:.3f}~cm$, $\sigma={bestfit_params.sigma:.3f}~cm$")
else:
    ax1.plot(z, bestfit_params.curve, color="black", linewidth = lineWidth, label=fr"No target")# $R_0={bestfit_params.R0:.3f}~cm$, $\sigma={bestfit_params.sigma:.3f}~cm$")
    #ax1.plot(z, depth_dose_distribution_ionization(z, bestfit_params.Phi0, bestfit_params.R0, bestfit_params.sigma, bestfit_params.epsilon), linewidth = 3, color="red")
    #ax1.plot(z, depth_dose_distribution_nonelastic(z, bestfit_params.Phi0, bestfit_params.R0, bestfit_params.sigma, bestfit_params.epsilon), linewidth = 3, color="blue")

    hist = ROOT.TH1D("h_fit", "heteroFit", len(z), z[0], z[-1])
    for i, value in enumerate(bestfit_params.curve):
        hist.SetBinContent(i+1, value)
    
    t=0
    sigmat=0
    pmod=0
    
target_fit_curves = []
target_fit_R0 = []
target_fit = []
pmma_thick = []

if(bhetero):
    spline_func = interp1d(depthTarget, unfoldedDoseTarget, kind='cubic')
    z_spline    = np.linspace(min(depthTarget), max(depthTarget), round((max(depthTarget)-min(depthTarget)) / resolution ))
    quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))
    
    R0Target = quantities['R80D']
    Phi0Target = max(unfoldedDoseTarget)/20
    sigmaMonoTarget = (beta*R0Target**0.935)
    sigmaTarget     = np.sqrt(sigmaMonoTarget**2+sigmaE0**2*a_h2o**2*p_h2o**2*beamEnergy**(2*p_h2o-2))

    convParams = bortfeld_fit(depthTarget, unfoldedDoseTarget, Phi0Target, R0Target, sigmaTarget, epsilon, None)

    print()
    print("Heteogeneous Target")
    print(f"Phi0: {convParams.Phi0} \pm {convParams.stddev[0]}")
    print(f"R0: {convParams.R0} \pm {convParams.stddev[1]}")
    print(f"sigma: {convParams.sigma} \pm {convParams.stddev[2]}")
    print(f"epsilon: {convParams.epsilon} \pm {convParams.stddev[3]}")

    params_list = []
    for i in range(nbOfFitsHetero):
        if i%100 == 0:
            print(f"Fit: {i}/{nbOfFitsHetero}")
        x_with_err = []
        y_with_err = []
        for index, mean in enumerate(depthTarget):
            x_with_err.append(gaussian_with_cutoff(mean, depthErrTarget[index]))
        for index, mean in enumerate(unfoldedDoseTarget):
            y_with_err.append(gaussian_with_cutoff(mean, unfoldedDoseErrTarget[index]))
        #params = bortfeld_fit_hetero(x_with_err, y_with_err, bestfit_params.Phi0, R0, sigma, bestfit_params.epsilon, weights)
        params = bortfeld_fit(x_with_err, y_with_err, Phi0Target, R0Target, sigmaTarget, epsilon, weights)

        if(params.Phi0 == 0):
            continue
        else:
            params_list.append(params)
            params_list[-1].curve = depth_dose_distribution(z, params.Phi0, params.R0, params.sigma, params.epsilon)
    
    convParams.curve = depth_dose_distribution(z, convParams.Phi0, convParams.R0, convParams.sigma, convParams.epsilon)
    
    hist = ROOT.TH1D("h_fit", "heteroFit", len(z), z[0], z[-1])
    for i, value in enumerate(convParams.curve):
        hist.SetBinContent(i+1, value)

    t = bestfit_params.R0-convParams.R0
    o_t = np.sqrt(bestfit_params.stddev[1]**2+convParams.stddev[1]**2)
    
    sigmat = np.sqrt(convParams.sigma**2-bestfit_params.sigma**2)
    
    o_sigmat = np.sqrt(convParams.sigma**2/(convParams.sigma**2-bestfit_params.sigma**2)*convParams.stddev[2]**2+bestfit_params.sigma**2/(convParams.sigma**2-bestfit_params.sigma**2)*bestfit_params.stddev[2]**2)
    pmod = sigmat**2/t*10000
    sigma_pmod = np.sqrt((2*sigmat/t**2*o_sigmat)**2+(sigmat**2/t**2*o_t)**2)*10000
    
    # ax1.plot(z, convParams.curve, color=targetColorMap[0], linewidth = lineWidth, label=fr"Hetero. target fit paramters: $\frac{{\Phi_0}}{{N_0}}={convParams.Phi0:.3f}~\frac{{1}}{{cm^2}}$," "\n" rf"$R_0={convParams.R0:.3f}~cm$, $\sigma={convParams.sigma:.3f}~cm$, " "\n" fr"$t= {t:.3f}" "\pm" fr"{o_t:.3f}~cm, \sigma={sigmat:.3f} \pm {o_sigmat:.3f}~cm,~P_{{mod}}={pmod:.3f} \pm {sigma_pmod:.3f}$" rf"$~\mu m$")
    ax1.plot(z, convParams.curve, color=targetColorMap[0], linewidth = lineWidth, label=fr"LN300 " rf"$R_0={convParams.R0:.3f}~cm$, $\sigma={convParams.sigma:.3f}~cm$, " "\n" fr"$t= {t:.3f}" "\pm" fr"{o_t:.3f}~cm, \sigma={sigmat:.3f} \pm {o_sigmat:.3f}~cm,~P_{{mod}}={pmod:.3f} \pm {sigma_pmod:.3f}$" rf"$~\mu m$")
elif(targetSelect == 1):
    for idx, thickness  in enumerate(targetThicknesses):
        spline_func = interp1d(depthTarget[idx], unfoldedDoseTarget[idx], kind='cubic')
        z_spline    = np.linspace(min(depthTarget[idx]), max(depthTarget[idx]), round((max(depthTarget[idx])-min(depthTarget[idx])) / resolution ))
        quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))

        R0Target = quantities['R80D']
        Phi0Target = max(unfoldedDoseTarget[idx])/20
        sigmaMonoTarget = (beta*R0Target**0.935)
        sigmaTarget     = np.sqrt(sigmaMonoTarget**2+sigmaE0**2*a_h2o**2*p_h2o**2*beamEnergy**(2*p_h2o-2))

        sigmaMono = (beta*R0**0.935)
        sigma     = np.sqrt(sigmaMono**2+sigmaE0**2*a_h2o**2*p_h2o**2*beamEnergy**(2*p_h2o-2))

        convParams = bortfeld_fit(depthTarget[idx], unfoldedDoseTarget[idx], Phi0Target, R0Target, sigmaTarget, epsilon, None)

        print(f"{thickness} mm PMMA Target")
        print(f"Phi0: {convParams.Phi0} \pm {convParams.stddev[0]}")
        print(f"R0: {convParams.R0} \pm {convParams.stddev[1]}")
        print(f"sigma: {convParams.sigma} \pm {convParams.stddev[2]}")
        print(f"epsilon: {convParams.epsilon} \pm {convParams.stddev[3]}")
        print()
        params_list = []
        for i in range(nbOfFitsHetero):
            if i%100 == 0:
                print(f"Fit: {i}/{nbOfFitsHetero}")
            x_with_err = []
            y_with_err = []
            for index, mean in enumerate(depthTarget[idx]):
                x_with_err.append(gaussian_with_cutoff(mean, depthErrTarget[idx][index]))
            for index, mean in enumerate(unfoldedDoseTarget):
                y_with_err.append(gaussian_with_cutoff(mean, unfoldedDoseErrTarget[idx][index]))
            #params = bortfeld_fit_hetero(x_with_err, y_with_err, bestfit_params.Phi0, R0, sigma, bestfit_params.epsilon, weights)
            params = bortfeld_fit(x_with_err, y_with_err, bestfit_params.Phi0, R0, sigma, bestfit_params.epsilon, weights)

            if(params.Phi0 == 0):
                continue
            else:
                params_list.append(params)
                params_list[-1].curve = depth_dose_distribution(z, params.Phi0, params.R0, params.sigma, params.epsilon)

        convParams.curve = depth_dose_distribution(z, convParams.Phi0, convParams.R0, convParams.sigma, convParams.epsilon)

        hist = ROOT.TH1D(f"h_fit{thickness}", f"homoFit{thickness}", len(z), z[0], z[-1])
        for i, value in enumerate(convParams.curve):
            hist.SetBinContent(i+1, value)

        t = bestfit_params.R0-convParams.R0
        o_t = np.sqrt(bestfit_params.stddev[1]**2+convParams.stddev[1]**2)
        pmma_thick.append([t, o_t])
        # sigmat = np.sqrt(convParams.sigma**2 - bestfit_params.sigma**2)

        # o_sigmat = np.sqrt(convParams.sigma**2/(convParams.sigma**2-bestfit_params.sigma**2)*convParams.stddev[2]**2+bestfit_params.sigma**2/(convParams.sigma**2-bestfit_params.sigma**2)*bestfit_params.stddev[2]**2)
        # pmod = sigmat**2/t*10000
        # sigma_pmod = np.sqrt((2*sigmat/t**2*o_sigmat)**2+(sigmat**2/t**2*o_t)**2)*10000

        # ax1.plot(z, convParams.curve, color="orange", linewidth = lineWidth, label=fr"Hetero. target fit paramters: $\frac{{\Phi_0}}{{N_0}}={convParams.Phi0:.3f}~\frac{{1}}{{cm^2}}$," "\n" rf"$R_0={convParams.R0:.3f}~cm$, $\sigma={convParams.sigma:.3f}~cm$, " "\n" fr"$t= {t:.3f}" "\pm" fr"{o_t:.3f}~cm, \sigma={sigmat:.3f} \pm {o_sigmat:.3f}~cm,~P_{{mod}}={pmod:.3f} \pm {sigma_pmod:.3f}$" rf"$~\mu m$")
        ax1.plot(z, convParams.curve, color=targetColorMap[idx], linewidth = lineWidth, label=fr"{thickness} mm PMMA")# " rf"$R_0={convParams.R0:.3f}~cm$, $\sigma={convParams.sigma:.3f}~cm$, " "\n" fr"$t= {t:.3f}" "\pm" fr"{o_t:.3f}~cm$")
        target_fit_curves.append(convParams.curve)
        target_fit_R0.append(convParams.R0)
        target_fit.append(convParams)
        
end_time = time.time()

elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")

ax1.grid(True)
ax1.set_xlabel('Water Equivalent Depth / cm')
ax1.set_ylabel('Dose  / μG')
ax1.legend(loc='upper left',  fancybox=False, edgecolor='black')
fig.tight_layout()

if(targetSelect > 0):
    plt.savefig(f"../data/{dataset}/{targetFile}/output/pdf/braggfit.pdf", format='pdf', bbox_inches='tight')
    plt.savefig(f"../data/{dataset}/{targetFile}/output/pdf/braggfit.svg", format='svg', bbox_inches='tight')
    plt.savefig(f"../data/{dataset}/{targetFile}/output/pdf/braggfit.jpg", format='jpg', bbox_inches='tight', dpi=600)
else:    
    print(f"Save to PDF: ../data/{dataset}/{file}/output/pdf/braggfit.pdf")
    plt.savefig(f"../data/{dataset}/{file}/output/pdf/braggfit.pdf", format='pdf', bbox_inches='tight')
    plt.savefig(f"../data/{dataset}/{file}/output/pdf/braggfit.svg", format='svg', bbox_inches='tight')
    plt.savefig(f"../data/{dataset}/{file}/output/pdf/braggfit.png", format='png', bbox_inches='tight', dpi=600)

zoom_half_width = 5  # cm (adjust to taste: 1–2 cm is typical)
zmin_zoom = 24 - 3
zmax_zoom = 25 + 1

# ============================================================
# Zoomed-in Bragg peak plot (targetSelect == 1)
# ============================================================

if targetSelect == 1:
    plt.rcParams.update({'font.size': 28})
    fig_zoom, axz = plt.subplots(figsize=(16, 12))

    # axz.set_title("Bragg peak region (zoomed)")
    axz.set_xlim(zmin_zoom, zmax_zoom)

    # --- No target: FULL error bars ---
    axz.errorbar(
        depth,
        unfoldedDose,
        unfoldedDoseErr,
        depthErr,
        fmt='s',
        markersize=3,
        capsize=capSize,
        elinewidth=lineWidth,
        color='black',
    )

    # --- No target fit ---
    axz.plot(
        z,
        bestfit_params.curve,
        color='black',
        linewidth=2,
        label=rf"No target"# ($R_0={bestfit_params.R0:.3f}$ cm)"
    )

    print(f"${bestfit_params.R0:.3f} \\pm {bestfit_params.stddev[1]:.3f}$ & ${bestfit_params.sigma:.3f} \\pm {bestfit_params.stddev[2]:.3f}$")

    # --- Homogeneous targets: data + fits ---
    for idx, thickness in enumerate(targetThicknesses):
        # data with error bars
        pmma_thick
        print(f"${target_fit[idx].R0:.3f} \\pm {target_fit[idx].stddev[1]:.3f}$ & ${target_fit[idx].sigma:.3f} \\pm {target_fit[idx].stddev[2]:.3f}$ & ${pmma_thick[idx][0]:.3f} \\pm {pmma_thick[idx][1]:.3f}$")

        
        axz.errorbar(
            depthTarget[idx],
            unfoldedDoseTarget[idx],
            unfoldedDoseErrTarget[idx],
            depthErrTarget[idx],
            fmt='o',
            markersize=3,
            capsize=capSize,
            elinewidth=lineWidth,
            color=targetColorMap[idx],
            alpha=0.9
        )

        # fitted curve
        axz.plot(
            z,
            target_fit_curves[idx],
            color=targetColorMap[idx],
            linewidth=2,
            linestyle='-',
            label=rf"{thickness} mm PMMA"# ($R_0={target_fit_R0[idx]:.3f}$ cm)"
        )

    axz.set_xlabel("Water Equivalent Depth / cm")
    axz.set_ylabel("Dose / μG")

    axz.grid(True)
    axz.legend(
        loc="upper left",
        fancybox=False,
        edgecolor="black"
    )

    fig_zoom.tight_layout()

    plt.savefig(
        f"../data/{dataset}/{targetFile}/output/pdf/braggfit_zoom.pdf",
        format="pdf",
        bbox_inches="tight"
    )
    
if plotEnable:
    plt.show()
if nosave: 
    exit()

Phi0_branch = np.zeros(1, dtype='float32')
R0_branch = np.zeros(1, dtype='float32')
sigma_branch = np.zeros(1, dtype='float32')
epsilon_branch = np.zeros(1, dtype='float32')

Pmod_branch = np.zeros(1, dtype='float32')
t_branch = np.zeros(1, dtype='float32')
sigmat_branch = np.zeros(1, dtype='float32')

curve_branch = np.zeros(len(z), dtype='float32')

if(bhetero):
    fit_file1 = ROOT.TFile(f"../data/{dataset}/{targetFile}/output/{targetFile}_fit.root", "RECREATE")
    vtree = ROOT.TTree("vtree", "Tree holding fit parameters")
    Phi0_branch[0] = convParams.Phi0
    R0_branch[0] = convParams.R0
    sigma_branch[0] = convParams.sigma
    epsilon_branch[0] = convParams.epsilon

    Pmod_branch[0] = pmod
    t_branch[0] = t
    sigmat_branch[0] = sigmat

    curve_branch[:] = convParams.curve
    hist.Write()
else:
    fit_file1 = ROOT.TFile(f"../data/{dataset}/{file}/output/{file}_fit.root", "RECREATE")
    vtree = ROOT.TTree("vtree", "Tree holding fit parameters")
    
    Phi0_branch[0] = bestfit_params.Phi0
    R0_branch[0] = bestfit_params.R0
    sigma_branch[0] = bestfit_params.sigma
    epsilon_branch[0] = bestfit_params.epsilon

    Pmod_branch[0] = pmod
    t_branch[0] = t
    sigmat_branch[0] = sigmat

    curve_branch[:] = bestfit_params.curve
    hist.Write()

vtree.Branch("Phi0", Phi0_branch, "Phi0/F")
vtree.Branch("R0", R0_branch, "R0/F")
vtree.Branch("sigma", sigma_branch, "sigma/F")
vtree.Branch("epsilon", epsilon_branch, "epsilon/F")
vtree.Branch("pmod", Pmod_branch, "pmod/F")
vtree.Branch("sigmat", sigmat_branch, "sigmat/F")
vtree.Branch("t", t_branch, "t/F")
vtree.Branch("curve", curve_branch, f"curve[{len(z)}]/F")

vtree.Fill()

for param in params_list:
    Phi0_branch[0] = param.Phi0
    R0_branch[0] = param.R0
    sigma_branch[0] = param.sigma
    epsilon_branch[0] = param.epsilon
    curve_branch[:] = param.curve
    vtree.Fill() #change

if fit_file1.IsOpen():
    print("Fit file is successfully written.")
else:
    print("Failed to write the fit file.")
fit_file1.Write()
fit_file1.Close()
