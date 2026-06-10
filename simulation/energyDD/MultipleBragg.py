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


def convolution_fit(Xconv, Yconv, params):
    popt, pcov = curve_fit(lambda x, amp, mean, stddev: depth_dose_convolved(x, params.Phi0, params.R0, params.sigma, params.epsilon, amp, mean, stddev), xdata= Xconv, ydata=Yconv, p0 = [1, 2, 0.2], bounds=((0.9, 0, 0), (1.1, 10, 2)))
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
    beta =  0.012
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

def bortfeld_fit(x, y, Phi0, R0, Sigma, epsilon):    
    params = fit_params()
    popt, pcov = curve_fit(lambda z, Phi0, R0, Sigma, epsilon: depth_dose_distribution(z, Phi0, R0, Sigma, epsilon), x, y, p0=[Phi0, R0, Sigma*10, 0.1*epsilon], bounds=((Phi0*0.25, R0 - 3*Sigma, 0.001*Sigma, 0), (Phi0*1.5, R0 + 3.5*Sigma, 10*Sigma, 2)), maxfev=int(1e8))
    params.curve = depth_dose_distribution(z, *popt)
    params.Phi0 = popt[0] 
    params.R0 = popt[1] 
    params.sigma = popt[2]
    params.epsilon = popt[3]
    params.stddev = np.sqrt(np.diag(pcov))
    return params

def bortfeld_fit_hetero(x, y, Phi0, R0, sigma, epsilon):    
    
    params = fit_params()
    #popt, _ = curve_fit(lambda z, R0, sigma: depth_dose_distribution(z, Phi0, R0, sigma, epsilon), x, y, p0=[R0, sigma], bounds=((R0 - 3*sigma, 0.1*sigma), (R0 + 3.5*sigma, 3*sigma)), maxfev=int(1e8))
    popt, pcov = curve_fit(lambda z, Phi0, R0, sigma: depth_dose_distribution(z, Phi0, R0, sigma, epsilon), x, y, p0=[Phi0, R0, sigma], bounds=((Phi0*0.5, R0 - 3*sigma, 0.01*sigma), (Phi0*1.5, R0 + 3.5*sigma, 3*sigma)), maxfev=int(1e8))
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

def rangeEnergy(E, alpha, p):
    return alpha*E**p

def Energyrange(R, alpha, p):
    return (R/alpha)**(1/p)

def S_Energy(E, alpha, p):
    return 1/(alpha*p)*E**(1-p)

def S_Range(R, alpha, p):
    return 1/(alpha**(1/p)*p)*R**(1/p-1)

def stragglingWidth(TargetThickness, initialEnergy):
    TargetThickness = np.array(TargetThickness)
    print(TargetThickness)
    alpha_h2o = 2.585e-2
    alpha_prime_h2o = 0.087/10
    p_h2o = 1.738

    alpha_pbwo4 = 7.275e-3
    alpha_prime_pbwo4 = 0.537/10
    p_pbwo4 = 1.77

    c1 = (alpha_prime_h2o*p_h2o**3*alpha_h2o**(2/p_h2o))/(3*p_h2o-2)
    c2 = (alpha_prime_pbwo4*p_pbwo4**3*alpha_pbwo4**(2/p_pbwo4))/(3*p_pbwo4-2)
    
    print(f"C1: {c1}, C2: {c2}")
    print(f"Sqrt C1: {np.sqrt(c1)}, C2: {np.sqrt(c2)}")
    R_h2o = rangeEnergy(initialEnergy, alpha_h2o, p_h2o)
    R_pbwo4 = rangeEnergy(initialEnergy, alpha_pbwo4, p_pbwo4)

    print(f"Ranges: R_h2o = {R_h2o} cm, R_pbwo4 = {R_pbwo4} cm")

    SToPWO = S_Energy(initialEnergy, alpha_h2o, p_h2o)/S_Energy(initialEnergy, alpha_pbwo4, p_pbwo4)
    SToH2O = S_Energy(initialEnergy, alpha_pbwo4, p_pbwo4)/S_Energy(initialEnergy, alpha_h2o, p_h2o)

    print(f"Stopping power ratios: SToPWO = {SToPWO}, SToH2O = {SToH2O}")

    # h2oVar = c1*((R_h2o)**(3-2/p_h2o)-(R_h2o-TargetThickness)**(3-2/p_h2o))
    # pbwo4Var = c2*((R_pbwo4-TargetThickness*SToPWO)**(3-2/p_pbwo4))*SToH2O

    onlyh2oVar = c1*((TargetThickness)**(3-2/p_h2o))
    onlypbwo4Var = c2*((TargetThickness*SToPWO)**(3-2/p_pbwo4))*SToH2O*SToH2O
    
    h2oVar = c1*((TargetThickness)**(3-2/p_h2o))
    pbwo4Var = c2*((R_pbwo4)**(3-2/p_pbwo4)-(TargetThickness*SToPWO)**(3-2/p_pbwo4))*SToH2O**2
    
    fullvarh2o = c1*((R_h2o)**(3-2/p_h2o))
    fullvarpwo = c2*((R_pbwo4)**(3-2/p_pbwo4))*SToH2O**2
    print(f"Fullsigma: {np.sqrt(fullvarh2o)/10}, PWO {np.sqrt(fullvarpwo)/10}")
    print(f"Variance contributions: h2oVar = {h2oVar}, pbwo4Var = {pbwo4Var*SToH2O}")

    sigma2 = h2oVar + pbwo4Var

    sigma = np.sqrt(np.array(sigma2))*1/10

    print(sigma)
    return sigma, np.sqrt(onlyh2oVar)*1/10, np.sqrt(onlypbwo4Var)*1/10

def straggling(z):

    alpha_h2o = 2.585e-2
    alpha_prime_h2o = 0.087/10
    p_h2o = 1.738

    alpha_pbwo4 = 7.275e-3
    alpha_prime_pbwo4 = 0.537/10
    p_pbwo4 = 1.690

    R_h2o = rangeEnergy(initialEnergy, alpha_h2o, p_h2o)
    R_pbwo4 = rangeEnergy(initialEnergy, alpha_pbwo4, p_pbwo4)

    SToPWO = S_Energy(initialEnergy, alpha_h2o, p_h2o)/S_Energy(initialEnergy, alpha_pbwo4, p_pbwo4)
    SToH2O = S_Energy(initialEnergy, alpha_pbwo4, p_pbwo4)/S_Energy(initialEnergy, alpha_h2o, p_h2o)

    c1 = (alpha_prime_h2o*(1/(p_h2o*alpha_h2o**(1/p_h2o))*(R_h2o-z)**(1/p_h2o-1))**(-2))
    c2 = (alpha_prime_pbwo4*(1/(p_pbwo4*alpha_pbwo4**(1/p_pbwo4))*(R_pbwo4-z*SToPWO)**(1/p_pbwo4-1))**(-2))*SToH2O
    
    print(f"C1: {c1}, C2: {c2}")
    print(f"Sqrt C1: {np.sqrt(c1)}, C2: {np.sqrt(c2)}")


    print(f"Ranges: R_h2o = {R_h2o} cm, R_pbwo4 = {R_pbwo4} cm")



    return c1, c2

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

datasets = ["data"]
in_data = ["notarget", "homotarget", "heterotarget"]
in_title = ["without a target", "with the homogeneous target", "with the heterogeneous target"]

dataset = datasets[0]
file = in_data[0]
targetFile = in_data[targetSelect]



lineWidth = 2
capSize = 3


sigmasPWO = []
sigmasH2O = []

materials = ["pbwo4", "h2o"]

for material in materials:
    depth = []
    depthTarget = []
    depthTargetTemp = []
    depthErr = []
    depthTargetErr = []
    depthTargetErrTemp = []
    dose = []
    doseTarget = []
    doseTargetTemp = []
    doseErr = []
    doseTargetErr = []
    doseTargetErrTemp = []
    unfoldedDose = []
    unfoldedDoseTarget = []
    unfoldedDoseTargetTemp = []
    unfoldedDoseErr = []
    unfoldedDoseTargetErr = []
    unfoldedDoseTargetErrTemp = []

    amplitudeTarget = []
    meanTarget = []
    sigmaTarget = []
    mean_errTarget = []
    sigma_errTarget = []

    unfoldedEntries = []
    unfoldedEntriesTarget = []
    unfoldedEntriesTargetTemp = []
    dataFile = np.load(f"{dataset}/{file}/input/depthdose_{material}_{nLayers}.npz")
    boolWET = True
    if material == "h2o":
        boolWET = True
    depth = dataFile["depth"]
    depthErr = dataFile["depth_err"]
    if boolWET:
        depth = dataFile["depthWET"]
        depthErr = dataFile["depthWET_err"]
    dose = dataFile["dose"] 
    doseErr = dataFile["dose_err"] 
    unfoldedEntries = dataFile["dose"][0]
    amplitude = dataFile["amplitude"]
    mean = dataFile["mean"]
    sigma = dataFile["sigma"]
    mean_err = dataFile["mean_err"]
    sigma_err = dataFile["sigma_err"]

    targetColorMap = ["#1f77b4", "#4e79a7", "#76b7b2", "#bab0ac", "#f28e2b", "#e15759", "#9c755f"]

    targetThicknesses = []
    for thick in range(0, 315, 5):
        targetThicknesses.append(thick)
    

    initialEnergy = Energyrange(31, a_h2o, p_h2o)
    print(f"Initial energy calculated from R0: {initialEnergy} MeV")
    sigma, h2ovar, pwovar = stragglingWidth([0]+targetThicknesses, initialEnergy)
    
    c1, c2 = straggling(np.array([0]+targetThicknesses))

    plt.figure(figsize=(12, 8))
    plt.plot([0]+targetThicknesses, c1, label="c1", marker='o')
    plt.plot([0]+targetThicknesses, c2, label="c2", marker='s')
    plt.legend()
    plt.grid()
    plt.show()
    plt.close()

    plt.figure(figsize=(12, 8))
    plt.plot([0]+targetThicknesses, sigma, label="sigma", marker='o')
    plt.plot([0]+targetThicknesses, h2ovar, label="H2Ovar", marker='s')
    plt.plot([0]+targetThicknesses, pwovar, label="pwovar", marker='^')
    plt.legend()
    plt.grid()
    plt.show()
    plt.close()
    
    # exit()

    targetThicknesses = []
    for thick in range(10, 200, 10):
        targetThicknesses.append(thick)

    for thickness in targetThicknesses:
        dataFile = np.load(f"{dataset}/{targetFile}/input/depthdose_{material}_{nLayers}_{thickness}.npz")
        if boolWET:
            depthTarget.append(dataFile["depthWET"] + thickness/10)
            depthTargetErr.append(dataFile["depthWET_err"])
        else:
            depthTarget.append(dataFile["depth"])
            depthTargetErr.append(dataFile["depth_err"])
        doseTarget.append(dataFile["dose"]) 
        doseTargetErr.append(dataFile["dose_err"])
        unfoldedEntriesTarget.append(dataFile["dose"][0])
        amplitudeTarget.append(dataFile["amplitude"])
        meanTarget.append(dataFile["mean"])
        sigmaTarget.append(dataFile["sigma"])
        mean_errTarget.append(dataFile["mean_err"])
        sigma_errTarget.append(dataFile["sigma_err"])


        # alpha = 0.02585
        # p = 1.738
        # t = mean-meanTarget
        # sigmat = np.sqrt(sigmaTarget**2-sigma**2)
        # pmod = sigmat**2/t
        # print("─────────────────────────────────────────────────────────")
        # print(f"Energy based pmod calculation")
        # print(f"t = {t} cm")
        # print(f"sigmat = {sigmat} cm")
        # print(f"pmod = {pmod} um")

        # t_E = alpha*(mean**p- meanTarget**p)
        # print(f"t (from Energy): {t_E} cm")
        # print("─────────────────────────────────────────────────────────")

    plt.rcParams.update({'font.size': 32})
    fig, ax1 = plt.subplots(figsize=(16, 12))

    start_time = time.time()

    z = np.linspace(0, 40, 4001)
    if material == "h2o":
        z = np.linspace(0, 40, 4001)
    else:
        z = np.linspace(0, 40, 4001)


    #─────────────────────────────────────────────────────────
    # Initial Guess calculation
    #─────────────────────────────────────────────────────────
    epsilon = 0.001*beamEnergy
    Phi0 = max(dose)/20
    beta = 0.012


    resolution = 0.01*np.min(np.diff(z))
    spline_func = interp1d(depth, dose, kind='cubic')
    z_spline    = np.linspace(min(depth), max(depth), round((max(depth)-min(depth)) / resolution ))
    quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))

    R0 = quantities['R80D']
    print("Expected range from R80D fit: ", range_energy_relationship(235, a_pwo, p_pwo))
    sigmaMono = (0.012*31**0.935)
    print("Expected sigma from monoenergetic beam: ", sigmaMono)
    # exit()
    sigmaE0   = 0.01*beamEnergy
    sigma     = np.sqrt(sigmaMono**2+sigmaE0**2*a_h2o**2*p_h2o**2*beamEnergy**(2*p_h2o-2))
    #─────────────────────────────────────────────────────────


    popt_notarget = bortfeld_fit(depth, dose, Phi0, R0, sigma, epsilon)

    ax1.errorbar(depth, dose, doseErr, depthErr, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color="#000000")#, label="No target data")
    for index, thickness in enumerate(targetThicknesses):
        ax1.errorbar(depthTarget[index], doseTarget[index], doseTargetErr[index], depthTargetErr[index], fmt='o', markersize=1, capsize=capSize, elinewidth=lineWidth, color=targetColorMap[index%len(targetColorMap)])#, label="Hetero. data") 

    print("─────────────────────────────────────────────────────────")
    print("Notarget fit parameters:")
    print("─────────────────────────────────────────────────────────")
    print(rf"Phi0: {popt_notarget.Phi0:.3f} \pm {popt_notarget.stddev[0]:.3f}")
    print(rf"R0: {popt_notarget.R0:.3f} \pm {popt_notarget.stddev[1]:.3f}")
    print(rf"sigma: {popt_notarget.sigma:.3f} \pm {popt_notarget.stddev[2]:.3f}")
    print(rf"epsilon: {popt_notarget.epsilon:.3f} \pm {popt_notarget.stddev[3]:.3f}")
    print("─────────────────────────────────────────────────────────")
    
    if material == "h2o":
        sigmasH2O.append(popt_notarget.sigma)
    else:
        sigmasPWO.append(popt_notarget.sigma)

    ax1.plot(z, popt_notarget.curve, color="black", linewidth = lineWidth, label=fr"No target")# rf"$R_0={popt_notarget.R0:.3f}~cm$, $\sigma={popt_notarget.sigma:.3f}~cm$")

    print("Heterotargets fit parameters:")
    print("─────────────────────────────────────────────────────────")
    for index, thickness in enumerate(targetThicknesses):
        #─────────────────────────────────────────────────────────
        # Initial Guess calculation
        #─────────────────────────────────────────────────────────
        spline_func = interp1d(depthTarget[index], doseTarget[index], kind='cubic')
        z_spline    = np.linspace(min(depthTarget[index]), max(depthTarget[index]), round((max(depthTarget[index])-min(depthTarget[index])) / resolution ))
        quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))

        R0Target = quantities['R80D']
        Phi0Target = max(doseTarget[index])/20
        sigmaMonoTarget = (beta*R0Target**0.935)
        sigmaTarget     = np.sqrt(sigmaMonoTarget**2+sigmaE0**2*a_h2o**2*p_h2o**2*beamEnergy**(2*p_h2o-2))

        popt_hetero = bortfeld_fit(depthTarget[index], doseTarget[index], Phi0Target, R0Target, sigmaTarget, epsilon)

        #─────────────────────────────────────────────────────────
        # Calculate Modulation power parameters and uncertainties
        #─────────────────────────────────────────────────────────
        t = popt_notarget.R0-popt_hetero.R0
        o_t = np.sqrt(popt_notarget.stddev[1]**2+popt_hetero.stddev[1]**2)

        sigmat = np.sqrt(popt_hetero.sigma**2-popt_notarget.sigma**2)

        o_sigmat = np.sqrt(popt_hetero.sigma**2/(popt_hetero.sigma**2-popt_notarget.sigma**2)*popt_hetero.stddev[2]**2+popt_notarget.sigma**2/(popt_hetero.sigma**2-popt_notarget.sigma**2)*popt_notarget.stddev[2]**2)
        pmod = sigmat**2/t*10000
        sigma_pmod = np.sqrt((2*sigmat/t**2*o_sigmat)**2+(sigmat**2/t**2*o_t)**2)*10000

        ax1.plot(z, popt_hetero.curve, color=targetColorMap[0], linewidth = lineWidth, label=fr"LN300 " rf"$P_{{mod}}={pmod:.3f} \pm {sigma_pmod:.3f}$" rf"$~\mu m$")
        print(f"Targetthickness {thickness} cm:")
        print("──────")
        print(rf"Phi0: {popt_hetero.Phi0:.3f} \pm {popt_hetero.stddev[0]:.3f}")
        print(rf"R0: {popt_hetero.R0:.3f} \pm {popt_hetero.stddev[1]:.3f}")
        print(rf"sigma: {popt_hetero.sigma:.3f} \pm {popt_hetero.stddev[2]:.3f}")
        print(rf"epsilon: {popt_hetero.epsilon:.3f} \pm {popt_hetero.stddev[3]:.3f}")
        print("──────")
        print(rf"t: {t:.3f} $\pm$ {o_t:.3f} cm")
        print(rf"sigmat: {sigmat:.3f} $\pm$ {o_sigmat:.3f} cm")
        print(rf"Pmod: {pmod:.3f} $\pm$ {sigma_pmod:.3f} um")
        print("─────────────────────────────────────────────────────────")
        if material == "h2o":
            sigmasH2O.append(popt_hetero.sigma)
        else:
            sigmasPWO.append(popt_hetero.sigma)
    end_time = time.time()

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    ax1.grid(True)
    ax1.set_xlabel('Water Equivalent Depth / cm')
    ax1.set_ylabel('rel. Dose  / μG')
    # ax1.legend(loc='upper left',  fancybox=False, edgecolor='black')
    fig.tight_layout()

    plt.savefig(f"{dataset}/{targetFile}/output/pdf/braggfit.pdf", format='pdf', bbox_inches='tight')
    # plt.show()
    plt.close()

initialEnergy = Energyrange(popt_notarget.R0, a_h2o, p_h2o)
print(f"Initial energy calculated from R0: {initialEnergy} MeV")
sigma, h2ovar, pwovar = stragglingWidth([0]+targetThicknesses, initialEnergy)

plt.figure(figsize=(16, 12))
plt.plot([0]+targetThicknesses, sigma, label="sigma", marker='o')
plt.plot([0]+targetThicknesses, h2ovar, label="H2Ovar", marker='s')
plt.plot([0]+targetThicknesses, pwovar, label="pwovar", marker='^')
plt.legend()
plt.grid()
plt.show()


sigmaDiff = np.array(sigmasH2O)/np.array(sigmasPWO)
plt.figure(figsize=(16, 12))
plt.plot([0]+targetThicknesses, sigmasH2O, label="H2O", marker='o')
plt.plot([0]+targetThicknesses, sigmasPWO, label="PWO", marker='s')
plt.plot([0]+targetThicknesses, sigmaDiff, label="Ratio", marker='^')
plt.plot([0]+targetThicknesses, sigma, label="Calculated Straggling Width", marker='x')
plt.xlabel("Target Thickness (mm)")
plt.ylabel("Sigma (cm)")
plt.title("Sigma vs Target Thickness")
plt.legend()
plt.savefig(f"{dataset}/{targetFile}/output/pdf/sigma_vs_thickness.pdf", format='pdf', bbox_inches='tight')
plt.show()
plt.close()