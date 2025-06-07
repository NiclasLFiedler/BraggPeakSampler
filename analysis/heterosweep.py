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

plt.style.use(['science','notebook','grid']) 

alpha  = 2.585e-2
alpha_o  = 6.826e-4
p = 1.738
p_o = 4.95e-3

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

def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-(x-mean)**2 / (2 * stddev**2)) / (np.sqrt(2 * np.pi) * stddev)

def GetRange(E):
    return alpha*E**p

def GetSigmaR(E, sigmaE):
    # print()
    # print(E)
    # print(sigmaE)
    # print(sigmaE*alpha*p*E**(p-1))
    # print(alpha_o*E**p)
    # print(alpha*E**p*np.log(E)*p_o)
    return np.sqrt((sigmaE*alpha*p*E**(p-1))**2)

def GetConv(E, sigma):
    return sigma*alpha*p*E**(p-1)

def range_energy_relationship(E, alpha, p):
    return alpha*E**p

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

def bortfeld_fit(x, y, Phi0, R0, epsilon, sigma, weights=None):    
    if weights is None:
        weights = np.ones_like(y)
    sigma_weights = 1 / weights
    
    params = fit_params()
    popt, pcov = curve_fit(lambda z, Phi0, R0, sigma, epsilon: depth_dose_distribution(z, Phi0, R0, sigma, epsilon), x, y, p0=[Phi0, R0, sigma*10, epsilon], bounds=((Phi0*0.5, R0 - 3*sigma, 0.1*sigma, 0), (Phi0*1.5, R0 + 3.5*sigma, 10*sigma, 1)), sigma=sigma_weights, maxfev=int(1e8))
    params.curve = depth_dose_distribution(x, *popt)
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
    params.curve = depth_dose_distribution(x, popt[0], popt[1], popt[2], epsilon)
    params.Phi0 = popt[0] 
    params.R0 = popt[1] 
    params.sigma = popt[2]
    params.epsilon = epsilon
    params.stddev = np.sqrt(np.diag(pcov))
    return params

colors = [
    "#1f77b4",  # muted blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # yellow-green
]


Meansfiles = []
Fitfiles = []
pmods = [100,200,300,400,500,600,700,800]
thicknesses = [50,100,150,200]

#pmods = [100, 200, 500]

#thicknesses = [200]

combination = []
combination.append([0,0,0])

for index, pmod in enumerate(pmods):
    for thickness in thicknesses:
        combination.append([pmod, thickness, index+1])
        
notargetR0 = 0
notargetR0sigma = 0

# Create a figure and axis for plotting
fig, ax = plt.subplots(figsize=(30, 13))

bPlotMeans = True
lineWidth = 2
capSize = 3

#output = "outputDratio"
output = "output"

notargetX = []
notargetY = []
z = np.linspace(0, 40, 4001)
f = []
f2 = []
popt = []
for comb in combination:
    meansfile = f"../data/heterosweep/{output}/{comb[0]}um_{comb[1]}mmMeans.root"
    fitfile = f"../data/heterosweep/{output}/{comb[0]}um_{comb[1]}mmFit.root"
    
    meansfileROOT = uproot.open(meansfile)
    fitfileROOT = ROOT.TFile(fitfile, "READ")
    
    meansTree = meansfileROOT["meantree"]

    y_data = meansTree["mean"].array().to_numpy()
    y_sigma = meansTree["error"].array().to_numpy()
    x_data = meansTree["x"].array().to_numpy()
    x_data = [x/10 for x in x_data]
    x_sigma = meansTree["x_sigma"].array().to_numpy()
    x_sigma = [x/10 for x in x_sigma]

    if bPlotMeans:
        ax.errorbar(x_data, y_data, y_sigma, x_sigma, fmt='s', markersize=1, capsize=capSize, elinewidth=lineWidth, color=colors[comb[2]])

    hist = fitfileROOT.Get("h_fit")
    fittree = fitfileROOT.Get("vtree")

    fittree.GetEntry(0)

    x = np.linspace(0,40,4001)
    beamEnergy = 220
    epsilon = 0.001*beamEnergy
    Phi0 = max(y_data)/20

    beta = 0.012

    resolution = 0.01*np.min(np.diff(x))

    spline_func = interp1d(x_data, y_data, kind='cubic')
    z_spline    = np.linspace(min(x_data), max(x_data), round((max(x_data)-min(x_data)) / resolution ))
    quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))


    R0 = quantities['R80D']
    if(comb[2] == 0):
        R0 = range_energy_relationship(beamEnergy, a_h2o, p_h2o)

    sigmaMono = (beta*R0**0.935)
    sigmaE0   = 0.01*beamEnergy
    sigma     = np.sqrt(sigmaMono**2+sigmaE0**2*a_h2o**2*p_h2o**2*beamEnergy**(2*p_h2o-2))

    if(comb[2] == 0):
        notargetX = x_data
        notargetY = y_data
        bestfit_params = bortfeld_fit(x_data, y_data, Phi0, R0, sigma, epsilon)
        t = 0
        o_t = 0
    
        sigmat = 0
        o_sigmat = 0

        pmod = 0
        sigma_pmod = 0
        params = bestfit_params
    else:
        convParams = bortfeld_fit(x_data, y_data, Phi0, R0, sigma, epsilon)
        params = convParams
        t = bestfit_params.R0-convParams.R0
        o_t = np.sqrt(bestfit_params.stddev[1]**2+convParams.stddev[1]**2)
    
        sigmat = np.sqrt(convParams.sigma**2 - bestfit_params.sigma**2)
        o_sigmat = np.sqrt(convParams.sigma**2/(convParams.sigma**2-bestfit_params.sigma**2)*convParams.stddev[2]**2+bestfit_params.sigma**2/(convParams.sigma**2-bestfit_params.sigma**2)*bestfit_params.stddev[2]**2)

        pmod = sigmat**2/t*10000
        sigma_pmod = np.sqrt((2*sigmat/t**2*o_sigmat)**2+(sigmat**2/t**2*o_t)**2)*10000

        f = interp1d(notargetX, notargetY, kind='linear', fill_value="extrapolate")
        f2 = interp1d(notargetX, notargetY, kind='cubic', fill_value="extrapolate")

        popt, pcov =  curve_fit(lambda x, amp, mean, stddev: right_sided_convolution(f, lambda x2: gaussian(x2, amp, mean, stddev), x), x_data, y_data, p0 = [1, 2, 0.2], bounds=((0.8, 0, 0), (1.3, 10, 2)))
        popt2, pcov2 =  curve_fit(lambda x, amp, mean, stddev: right_sided_convolution(f2, lambda x2: gaussian(x2, amp, mean, stddev), x), x_data, y_data, p0 = [*popt], bounds=((0.4, 0, 0), (1.3, 10, 2)))
        #print(popt)
        t_conv = popt[1]
        sigmat_conv = popt[2]
        pmod_conv = popt[2]**2/popt[1]*10**4
        
        t_conv2 = popt2[1]
        sigmat_conv2 = popt2[2]
        pmod_conv2 = popt2[2]**2/popt2[1]*10**4

    # Step 2: Plot using Matplotlib
    if comb[2] != 0:
        labeltext = f"{comb[0]} um, {comb[1]} mm, Fit: {params.R0:.3f} mm, {params.sigma:.3f} mm, {t:.3f} cm, {sigmat:.3f} cm, {pmod:.3f} um | Conv: {t_conv:.3f} cm, {sigmat_conv:.3f} cm, {pmod_conv:.3f} um, | Conv2: {t_conv2:.3f} cm, {sigmat_conv2:.3f} cm, {pmod_conv2:.3f} um ||  Diff: {(pmod/comb[0]-1)*100:.2f} % | {(pmod_conv/comb[0]-1)*100:.2f} % | {(pmod_conv2/comb[0]-1)*100:.2f} %"
        ax.plot(z, right_sided_convolution(f2, lambda x2: gaussian(x2, *popt), z), label='Right-sided convolved (Gaussian)')
    else:
        labeltext = f"{comb[0]} um, {comb[1]} mm, {params.R0:.3f} mm, {params.sigma:.3f} mm, {t:.3f} cm, {sigmat:.3f} cm, {pmod:.3f} um, Diff: {0} %"
    
    #ax.step(x, depth_dose_distribution(x, params.Phi0, params.R0, params.sigma, params.epsilon), where='mid', label=labeltext, linewidth=1, color=colors[comb[2]])
    
    print(labeltext)
    # texText = f"{comb[0]} & {comb[1]} & {popt[1]:.2f}" " $\\pm$ " f"{popt[2]:.2f} & {R0:.2f}" " $\\pm$ " f"{sigmaR0:.2f} & {Pmod:.2f} \\\\"
    # print(texText)

ax.plot(notargetX, f(notargetX))
#ax.set_yscale('log')
#ax.set_xlim(175, 225)
# Add labels and title
ax.set_xlabel('Energy / MeV')
ax.set_ylabel('Counts')
ax.set_title('Modulation power analysis from total energy distribution')

# Add a legend
#ax.legend(fontsize=12)
#fig.tight_layout()
# Show the plot
plt.show()
