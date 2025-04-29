import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pybragg import fitBP
from pybragg import bortfeld
import scienceplots
import uproot
from scipy.integrate import quad
import time
import pandas as pd
plt.style.use(['science','notebook','grid'])
import ROOT

class fit_params:
    def __init__(self, D100=0, R0=0, sigma=0, p=0, k=0, curve=[], perr = []) -> None:
        self.D100 = D100
        self.R0 = R0
        self.sigma = sigma
        self.p = p
        self.k = k
        self.curve = curve
        self.perr = perr

def gaussian_with_cutoff(mean, sigma, cutoff=2.5):
    while True:
        value = np.random.normal(mean, sigma)
        if mean-cutoff <= value <= mean+cutoff:
            return value

def bortfeld_fit(x, y):
    params = fit_params()
    try:
        bortfeld_fit_params = fitBP(x, y) # fit and get curve
    except:
        return params
    params.curve = bortfeld(z, *bortfeld_fit_params['bortfeld_fit_p'])
    params.perr = np.sqrt(np.diag(bortfeld_fit_params['bortfeld_fit_cov']))
    params.D100 = bortfeld_fit_params['bortfeld_fit_p'][0] 
    params.R0 = bortfeld_fit_params['bortfeld_fit_p'][1] 
    params.sigma = bortfeld_fit_params['bortfeld_fit_p'][2]
    params.p = bortfeld_fit_params['bortfeld_fit_p'][3]
    params.k = bortfeld_fit_params['bortfeld_fit_p'][4]
    return params

convolve_target = "beam2"
conv_file = uproot.open(f"MIT_05_2024/bp-p-{convolve_target}/output/{convolve_target}Means.root")
conv_tree = conv_file["meantree"]

conv_y = np.array(conv_tree["mean"].array())

conv_x = []
for value in np.array(conv_tree["ch"].array()):
    conv_x.append(0.25 + 0.8*value)

convolve_params = fitBP(conv_x, conv_y)
perr = np.sqrt(np.diag(convolve_params['bortfeld_fit_cov']))

print("Parameters:", convolve_params['bortfeld_fit_p'])
print("Parameter errors:", perr)

def gauss_function(x, mu, sigma):
     return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def convolve_integrand(x,sigma,t):
    return bortfeld(x, *convolve_params['bortfeld_fit_p'])*gauss_function(x, t, sigma)
    
def gaus_convolved_bortfeld(x0, sigma, t):
    gauss_kernel = gauss_function(x0, t, sigma)
    gauss_kernel /= np.sum(gauss_kernel)
    bortfeld_result = bortfeld(x0, *convolve_params['bortfeld_fit_p'])
    convolved_result = np.convolve(bortfeld_result, gauss_kernel, mode='same')[:len(x0)]
    return convolved_result

# Specify the file name
targets = ["beam2","homo-target","hetero-target"]
# file = uproot.open("MIT_05_2024/bp-p-beam1/output/histMeans.root")
# root_file = ROOT.TFile("curve_beam1.root", "RECREATE")

y_data = []
x_data = []
y_sigma = []
x_sigma = 0.5/np.sqrt(12)
file = []

plt.figure(figsize=(16, 9))
plt.title(f'Beamtime Energy dose of protons in PbWO4') 
z = np.linspace(0, 12, 1201) # depth in cm

for target in targets:  
    file = uproot.open(f"MIT_05_2024/bp-p-{target}/output/{target}Means.root")
    root_file = ROOT.TFile(f"MIT_05_2024/bp-p-{target}/curve_{target}.root", "RECREATE")

    tree = file["meantree"]
    
    y_data.append(np.array(tree["mean"].array()))
    y_sigma.append(np.array(tree["error"].array()))
    
    x_temp = []
    for value in np.array(tree["ch"].array()):
        x_temp.append(0.25 + 0.8*value)
    x_data.append(x_temp)
    
    

params = []
for index in range(3):
    params.append(bortfeld_fit(x_data[index], y_data[index]))

norm = params[0].curve.max()

for index in range(1,3):
    if(norm < params[index].curve.max()):
        norm = params[index].curve.max()

plt.xlabel('Depth / cm')
plt.ylabel('Normed Total Energy Dose')

#notarget

colorlist = ["blue","orange","green"]

t = (params[0].R0 - params[2].R0)*10
pmod = 0.256*(1-0.25)*(1/1)*1e3
sigma = np.sqrt(pmod/1e3*t/10)
for index in range(3):
    if(index == 2):
        plt.errorbar(x_data[index], y_data[index]/norm, y_sigma[index]/norm, x_sigma,fmt='o', color=colorlist[index], label=f"{targets[index]} data\n$R_0$={params[index].R0:.2e}~mm ± {params[index].perr[1]:.2e} mm\nt={t:.2e} mm\n"r"$P_{\text{mod}}"f"={pmod:.2e}~\mu m$\n$\sigma_t$ = {sigma:.2e} mm")
    else:
        plt.errorbar(x_data[index], y_data[index]/norm, y_sigma[index]/norm, x_sigma,fmt='o', color=colorlist[index], label=f"{targets[index]} data\n$R_0$={params[index].R0:.2e} mm ± {params[index].perr[1]:.2e} mm")
    plt.plot(z, params[index].curve/norm, color=colorlist[index])

# conv_p, conv_cov = curve_fit(gaus_convolved_bortfeld, x_data[2], y_data[2], maxfev=10000)
# print("conv_p")
# print(conv_p)
# conv_curve = gaus_convolved_bortfeld(z,*conv_p)
# #conv_curve = gaus_convolved_bortfeld(z, 0.0001, 4.2)
# print(conv_curve)
# plt.plot(z, conv_curve)#/conv_curve.max()


plt.legend(loc='upper left',  fancybox=False, edgecolor='black')
print("---------------------------------------")
print("        Heterogeneous analysis         ")
print("---------------------------------------")
print(f"No target R0: {(params[0].R0*10):.2e} mm")
print(f"Heterogeneous target R0: {(params[2].R0*10):.2e} cm")
print(f"Water equivalent thickness t: {t:.2e} mm")
# pmod = d * (1 - p_l) * rho_l/rho_h2o
print(f"Modulation Power Pmod: {pmod:.2e} mm")
print(f"Gaus distribution width Sigma: {sigma:.2e} mm")
print("---------------------------------------")
plt.savefig(f'MIT_05_2024/all_hist.pdf')
#plt.show()