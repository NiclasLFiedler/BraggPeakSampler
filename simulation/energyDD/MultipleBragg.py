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

dataFile = np.load(f"{dataset}/{file}/input/depthdose_{detectorType}_{nLayers}.npz")
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
for thick in range(10, 210, 10):
    targetThicknesses.append(thick)
 
for thickness in targetThicknesses:
    dataFile = np.load(f"{dataset}/{targetFile}/input/depthdose_{detectorType}_{nLayers}_{thickness}.npz")
    if boolWET:
        depthTarget.append(dataFile["depthWET"])
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
if detectorType == "h2o":
    z = np.linspace(0, crystalSize[2], 4001)
else:
    z = np.linspace(0, crystalSize[2]*5, 4001)


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
sigmaMono = (beta*R0**0.935)
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
        
end_time = time.time()

elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")

ax1.grid(True)
ax1.set_xlabel('Water Equivalent Depth / cm')
ax1.set_ylabel('rel. Dose  / μG')
# ax1.legend(loc='upper left',  fancybox=False, edgecolor='black')
fig.tight_layout()

plt.savefig(f"{dataset}/{targetFile}/output/pdf/braggfit.pdf", format='pdf', bbox_inches='tight')
plt.show()

Phi0_branch = np.zeros(1, dtype='float32')
R0_branch = np.zeros(1, dtype='float32')
sigma_branch = np.zeros(1, dtype='float32')
epsilon_branch = np.zeros(1, dtype='float32')

Pmod_branch = np.zeros(1, dtype='float32')
t_branch = np.zeros(1, dtype='float32')
sigmat_branch = np.zeros(1, dtype='float32')

curve_branch = np.zeros(len(z), dtype='float32')

if(bhetero):
    fit_file1 = ROOT.TFile(f"{dataset}/{targetFile}/output/{targetFile}_fit.root", "RECREATE")
    vtree = ROOT.TTree("vtree", "Tree holding fit parameters")
    Phi0_branch[0] = popt_hetero.Phi0
    R0_branch[0] = popt_hetero.R0
    sigma_branch[0] = popt_hetero.sigma
    epsilon_branch[0] = popt_hetero.epsilon

    Pmod_branch[0] = pmod
    t_branch[0] = t
    sigmat_branch[0] = sigmat

    curve_branch[:] = popt_hetero.curve

vtree.Branch("Phi0", Phi0_branch, "Phi0/F")
vtree.Branch("R0", R0_branch, "R0/F")
vtree.Branch("sigma", sigma_branch, "sigma/F")
vtree.Branch("epsilon", epsilon_branch, "epsilon/F")
vtree.Branch("pmod", Pmod_branch, "pmod/F")
vtree.Branch("sigmat", sigmat_branch, "sigmat/F")
vtree.Branch("t", t_branch, "t/F")
vtree.Branch("curve", curve_branch, f"curve[{len(z)}]/F")

vtree.Fill()


if fit_file1.IsOpen():
    print("Fit file is successfully written.")
else:
    print("Failed to write the fit file.")
fit_file1.Write()
fit_file1.Close()
