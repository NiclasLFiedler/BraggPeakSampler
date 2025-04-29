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

class fit_params:
    def __init__(self, Phi0=0, R0=0, sigma=0, epsilon=0, curve=[]) -> None:
        self.Phi0 = Phi0
        self.R0 = R0
        self.sigma = sigma
        self.epsilon = epsilon
        self.curve = curve

class fit_params_conv:
    def __init__(self, t=0, sigma=0, curve=[], cov = []) -> None:
        self.t = t
        self.sigma = sigma
        self.curve = curve
        self.cov = cov

# Define a Gaussian function
def gaussian(x, amp, mean, stddev):
    return np.exp(-(x-mean)**2 / (2 * stddev**2)) / (np.sqrt(2 * np.pi) * stddev)

def right_sided_convolution(f, g, z_values):
    """
    Compute the right-sided convolution (f * g)(z) = ∫ f(z') g(z'+z) dz'.

    Parameters:
        f (callable): Function f(z').
        g (callable): Function g(z'+z).
        z_values (array): Values of z at which to compute the convolution.

    Returns:
        array: Convolution results at each z in z_values.
    """
    def convolve_at_z(z):
        integrand = lambda z_prime: f(z_prime+z) * g(z_prime)
        result, _ = quad(integrand, 0, 40, limit=50)
        return result

    def convole_at_z2(z):
        dz = z[1]-z[0]
        convolution_result = np.zeros_like(z)
    
        for i, z_val in enumerate(z):
            shifted_f = lambda z: f(z + z_val)
            product = g(z) * shifted_f(z)          # Element-wise multiplication
            convolution_result[i] = simpson(product, dx=dz)  # Integrate over z using Simpson's rule

        return convolution_result

    convolve_func = np.vectorize(convolve_at_z)
    #return np.array([convole_at_z2(z) for z in z_values])
    return convole_at_z2(z_values)

def numerical_convolution(f_vals, g_vals, z_values, dx):
    conv_result = np.zeros_like(z_values)
    for i, z in enumerate(z_values):
        integrand = f_vals * np.roll(g_vals, int(z / dx))
        conv_result[i] = np.sum(integrand) * dx
    return conv_result

def array_to_callable(array, x_values):
    """Convert a NumPy array to a callable function using interpolation."""
    return interp1d(x_values, array, bounds_error=False, fill_value=0)

# Define a function for the Gaussian convolution
def depth_dose_convolved(x, Phi0, R0, sigma, epsilon, mean, stddev):
    x1=np.linspace(0,40,4001)
    indices = [np.abs(x1 - xi).argmin() for xi in x]

    reference_curve = depth_dose_distribution(x, Phi0, R0, sigma, epsilon)
    reference_curve = np.nan_to_num(reference_curve, nan=0, posinf=0, neginf=0)
    gaussian_kernel = gaussian(x, 1, mean, stddev)
    gaussian_kernel = np.nan_to_num(gaussian_kernel, nan=0, posinf=0, neginf=0)
    g_flipped = gaussian_kernel#[::-1]
    r_flipped = reference_curve[::-1]
    #convolved_curve = convolve(r_flipped, g_flipped/np.sum(g_flipped), mode='same')
    f_callable = array_to_callable(reference_curve, x)
    g_callable = array_to_callable(gaussian_kernel, x)
    convolved_curve = right_sided_convolution(f_callable, g_callable, x)
    #convolved_curve = ndimage.convolve(reference_curve, gaussian_kernel/np.sum(gaussian_kernel), mode='constant', cval=0.0)
    #convolved_curve = convolved_curve[(len(convolved_curve) - len(x)) // 2 : (len(convolved_curve) + len(x)) // 2]
    #return convolved_curve[indices]
    return convolved_curve

def convolution_fit(x_data,y_data, params, weights=None):
    
    # If weights are not provided, use uniform weights
    if weights is None:
        weights = np.ones_like(y_data)
    
    # Calculate sigma as the inverse of weights
    sigma = 1 / weights
    
    popt, pcov = curve_fit(lambda x, mean, stddev: depth_dose_convolved(x, params.Phi0, params.R0, params.sigma, params.epsilon, mean, stddev), x_data, y_data, p0 = [6.3, 0.30], sigma=sigma, bounds=((0, 0), (10, 1)), maxfev = 1000000)
    errors = np.sqrt(np.diag(pcov))
    return [popt, errors]

def cyl_gauss(a,x):
    "Calculate product of Gaussian and parabolic cylinder function"
    xarr=np.copy(x) # to make xarr a numpy array, even for scalar (0-D) arguments
    yarr = np.copy(xarr)
    branch1 = -12.0   # for large negative values of the argument we run into numerical problems, need to approximate result
    # branch1 = -1000.0
    branch2 = 20.0   # above branch2 function is sufficiently close to 0

    #indset1 = np.where(xarr<branch1)
    indset1 = (xarr<branch1) 
    x1 = xarr[indset1]
    y1 = np.sqrt(2*np.pi)/special.gamma(-a)*(-x1)**(-a-1)
    yarr[indset1] = y1

    #indset2 = np.where((xarr>=branch1) & (xarr<branch2))
    indset2 = ((xarr>=branch1) & (xarr<branch2))
    x2 = xarr[indset2]
    y2a = special.pbdv(a,x2)[0]     #special function yielding parabolic cylinder function, first array [0] is function itself
    y2b = np.exp(-x2*x2/4)
    y2 = y2a*y2b
    yarr[indset2] = y2

    #indset3 = np.where(xarr>=branch2)
    indset3 = (xarr>=branch2)
    yarr[indset3] = 0.0

    return yarr

def depth_dose_distribution(z, Phi0, R0, sigma, epsilon): 
    beta =  0
    gamma = 0.6
    p = 1.738
    a=2.585e-3

    return Phi0*sigma**(1/p)*special.gamma(1/p)/(np.sqrt(2*np.pi)*p*a**(1/p)*(1+beta*R0))*(sigma**(-1)*cyl_gauss(-1/p,(z-R0)/sigma)+(beta/p + gamma*beta + epsilon/R0)*cyl_gauss(-1/p-1,(z-R0)/sigma)) # eq26

def gaussian_with_cutoff(mean, sigma, cutoff=2.5):
    while True:
        value = np.random.normal(mean, sigma)
        if mean-cutoff <= value <= mean+cutoff:
            return value

def range_energy_relationship(E, alpha, p):
    return alpha*E**p

def bortfeld_fit(x, y, Phi0, R0, epsilon, sigma, weights=None):
    # If weights are not provided, use uniform weights
    if weights is None:
        weights = np.ones_like(y)
    
    # Calculate sigma as the inverse of weights
    sigma_weights = 1 / weights
    
    
    params = fit_params()
    popt, pcov = curve_fit(
        lambda z, Phi0, R0, sigma: depth_dose_distribution(z, Phi0, R0, sigma, 0),
        x, y,
        p0=[Phi0, R0, sigma,],  # Initial guesses for the four parameters
        bounds=((Phi0*0.5, R0 - 3*sigma, 0.1*sigma), (Phi0*1.5, R0 + 3.5*sigma, 3*sigma)),  # Bounds for the parameters
        sigma=sigma_weights,  # Weights for the fitting
        maxfev=int(1e8)  # Maximum function evaluations
    )
    params.curve = depth_dose_distribution(z, *popt,0)
    params.Phi0 = popt[0] 
    params.R0 = popt[1] 
    params.sigma = popt[2]
    params.epsilon = 0
    return params

def characterize_z_D_curve(z, D):
    """ Method that computes dose and range quantities from a given bragg curve
    
    Parameters
    -----------
    :param z: depth in phantom in cm
    :param D: dose at depth z
    
    Returns
    --------
    :returns: 
      - results (dict): Ranges to certain fractions of Dmax on distal (D) and proximal (P) side of peak. Also: FWHM, DFO(1090)/(2080)
    """
    
    results = {}

    # compute quantities
    D100_index = np.argmax(D)
    D100       = D[D100_index]
    R100       = z[D100_index]
    results.update({
        'D100': D100,
        'R100': R100
    })

    # split at peak index into proximal and distal part
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
    results.update({
        'R90P': R90P,
        'R90D': R90D,
        'R80P': R80P,
        'R80D': R80D,
        'R50P': R50P,
        'R50D': R50D,
        'R20D': R20D,
        'R10D': R10D
    })

    FWHM    = R50D  - R50P
    DFO2080 = R20D  - R80D
    DFO1090 = R10D  - R90D
    results.update({
        'FWHM': FWHM,
        'DFO2080': DFO2080,
        'DFO1090': DFO1090,
    })

    return results

# Specify the file name
beamtime = "MIT_05_2024"
target1 = "beam2" #"beam2" "homo-target"
target2 = "hetero-target"
hetero = True
nosave = True

file1 = uproot.open(f"{beamtime}/bp-p-{target1}/output/{target1}Means.root")

file2 = uproot.open(f"{beamtime}/bp-p-{target2}/output/{target2}Means.root")

#plt.style.use(['science','notebook','grid'])
plt.rcParams.update({'font.size': 20})
fig, ax1 = plt.subplots(figsize=(16, 8))
ax1.set_title(f'Fitted energy depth dose distribution without a target')
if(target1 == "homo-target"):
    print('Fitted energy depth dose distribution with the homogeneous target')
    ax1.set_title(f'Fitted energy depth dose distribution with the homogeneous target')
if(hetero):
    print('Fitted energy depth dose distribution with the heterogeneous target')
    ax1.set_title(f'Fitted energy depth dose distribution with the heterogeneous target')

tree1 = file1["meantree"]
tree2 = file2["meantree"]

z = np.linspace(0, 40, 4001) # depth in cm
z2 = np.linspace(0, 40, 401) # depth in cm
x_data = []
x_sigma = []

y_data1 = tree1["mean"].array().to_numpy()
y_data2 = tree2["mean"].array().to_numpy()

a_h2o=2.585e-3
a_pwo=7.275e-4
p_h2o=1.738
p_pwo=1.690
E = 220
wet_conv = (range_energy_relationship(E,a_h2o,p_h2o)/range_energy_relationship(E,a_pwo,p_pwo))
wet_conv2 = a_h2o/a_pwo*(p_h2o-p_pwo)*E**(p_h2o-p_pwo-1)
print(f"wet_conv: {wet_conv}")
#print(f"wet_conv2: {wet_conv2}")
for ch in tree1["ch"].array().to_numpy():
    x_data.append((0.25 + 0.5*ch)*wet_conv)
    x_sigma.append(np.sqrt((0.5/np.sqrt(12)*wet_conv)**2 + ((0.25 + 0.5*ch)*wet_conv2*0.005*E)**2))

y_sigma1 = tree1["error"].array().to_numpy()
y_sigma2 = tree2["error"].array().to_numpy()

start_time = time.time()

epsilon = 0.001*E
Phi0 = max(y_data1)/20
print(f"Phi0 {Phi0}")
beta = 0

resolution = 0.01*np.min(np.diff(z))

# fit spline with given precision to curve
spline_func = interp1d(x_data, y_data1, kind='cubic')
z_spline    = np.linspace(min(x_data), max(x_data), round((max(x_data)-min(x_data)) / resolution ))
quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))


R0 = quantities['R80D']
if target1 == "beam2":
    R0 = range_energy_relationship(E, a_h2o, p_h2o)

sigmaMono = (0.012*R0**0.935) # paper: Eq. (18)
sigmaE0   = 0.01*E # assumtion that Delta E will be small
sigma     = np.sqrt(sigmaMono**2+sigmaE0**2*a_h2o**2*p_h2o**2*E**(2*p_h2o-2))
print(sigma)
#sigma = np.sqrt((beta*R0**0.935)**2+(0.01*E*a_h2o*p_h2o)**2*E**(2*p_h2o-2))

top_indices = np.argsort(y_data1)[-3:]
weights = np.ones_like(y_data1)
weights[top_indices] = 1

nbOfFits = 0
params_list = []
for i in range(nbOfFits):
    if i%100 == 0:
        print(f"Fit: {i}/{nbOfFits}")
    x_with_err = []
    y_with_err = []
    for index, mean in enumerate(x_data):
        x_with_err.append(gaussian_with_cutoff(mean, x_sigma[index]))
    for index, mean in enumerate(y_data1):
        y_with_err.append(gaussian_with_cutoff(mean, y_sigma1[index]))
    params = bortfeld_fit(x_with_err, y_with_err, Phi0, R0, epsilon, sigma, weights)
    
    if(params.Phi0 == 0):
        continue
    else:
        params_list.append(params)
        params_list[-1].curve = depth_dose_distribution(z, params.Phi0, params.R0, params.sigma, params.epsilon)
        #plt.plot(z, params_list[-1].curve)
if(hetero):
    ax1.errorbar(x_data, y_data2, y_sigma2, x_sigma, fmt='o', markersize=2, capsize=3, elinewidth=1, color="black", label="Measured data points") #Convolution
else:
    ax1.errorbar(x_data, y_data1, y_sigma1, x_sigma, fmt='o', markersize=2, capsize=3, elinewidth=1, color="black", label="Measured data points") 
bestfit_params = bortfeld_fit(x_data, y_data1, Phi0, R0, epsilon, sigma, weights)
print(f"Phi0: {bestfit_params.Phi0}")
print(f"R0: {bestfit_params.R0}")
print(f"sigma: {bestfit_params.sigma}")
if(hetero):
    ax1.plot(z, bestfit_params.curve, color="grey", label=fr"Reference Curve: $\frac{{\Phi_0}}{{N_0}}={bestfit_params.Phi0:.3f}~\frac{{1}}{{cm^2}}$," "\n" rf"$R_0={bestfit_params.R0:.3f}~cm$, $\sigma={bestfit_params.sigma:.3f}~cm$")
else:
    ax1.plot(z, bestfit_params.curve, color="orange", label=fr"Fit: $\frac{{\Phi_0}}{{N_0}}={bestfit_params.Phi0:.3f}~\frac{{1}}{{cm^2}}$, $R_0={bestfit_params.R0:.3f}~cm$, $\sigma={bestfit_params.sigma:.3f}~cm$")

# Convolution 
top_indices = np.argsort(y_data2)[-3:]
weights = np.ones_like(y_data2)
weights[top_indices] = 1
[best_popt, best_error] = convolution_fit(x_data, y_data2, bestfit_params, weights)

bestfit_params_conv = fit_params_conv()
bestfit_params_conv.t = best_popt[0]
bestfit_params_conv.sigma = best_popt[1]
bestfit_params_conv.curve = depth_dose_convolved(z, bestfit_params.Phi0, bestfit_params.R0, bestfit_params.sigma, bestfit_params.epsilon, *best_popt)
bestfit_params_conv.cov = best_error

nbOfFits_conv = 0
params_list_conv = []
param_conv = fit_params_conv()
for i in range(nbOfFits_conv):
    if i%100 == 0:
        print(f"Fit: {i}/{nbOfFits_conv}")
    x_with_err = []
    y2_with_err = []
    for index, mean in enumerate(x_data):
        x_with_err.append(gaussian_with_cutoff(mean, x_sigma[index]))
    for index, mean in enumerate(y_data2):
        y2_with_err.append(gaussian_with_cutoff(mean, y_sigma2[index]))
    [popt_conv, pcov_conv] = convolution_fit(x_with_err, y2_with_err, bestfit_params, weights)
    param_conv.t = popt_conv[0]
    param_conv.sigma = popt_conv[1]
    param_conv.curve = depth_dose_convolved(z, bestfit_params.Phi0, bestfit_params.R0, bestfit_params.sigma, bestfit_params.epsilon, *popt_conv)
    params_list_conv.append(param_conv)

#print(best_popt)
if(hetero):
    #plt.plot(z2, depth_dose_convolved(z2, bestfit_params.Phi0, bestfit_params.R0, bestfit_params.sigma, bestfit_params.epsilon, 5, 0.3), color="orange")
    ax1.plot(z2, depth_dose_convolved(z2, bestfit_params.Phi0, bestfit_params.R0, bestfit_params.sigma, bestfit_params.epsilon, *best_popt), color="orange", label=fr"Convolved Fit: $t={best_popt[0]:.3f}~cm\pm{best_error[0]:.3f}~cm$," "\n" fr"$\sigma_t={best_popt[1]:.3f}~cm\pm{best_error[1]:.3f}~cm$," "\n" fr"$P_{{mod}}={(10000*best_popt[1])**2/(10000*best_popt[0]):.3f}~\mu m \pm{np.sqrt(((2*10000*best_popt[1])/(10000*best_popt[0])*10000*best_error[1])**2+(((10000*best_popt[1])**2)/((10000*best_popt[0])**2)*10000*best_error[0])**2):.3f}~\mu m$")

end_time = time.time()

elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")
#plt.legend()
ax1.grid(True)
ax1.set_xlabel('Depth / cm')
ax1.set_ylabel('Normalized Energy Dose / MeV')
fig.tight_layout()
ax1.legend(loc='upper left',  fancybox=False, edgecolor='black')
if(hetero):
    plt.savefig(f"{beamtime}/bp-p-{target2}/output/braggfit.pdf", format='pdf', bbox_inches='tight')
else:    
    plt.savefig(f"{beamtime}/bp-p-{target1}/output/braggfit.pdf", format='pdf', bbox_inches='tight')
if nosave:
    plt.show()
    exit()
# Create TTree
if(hetero):
    fit_file1 = ROOT.TFile(f"{beamtime}/bp-p-{target2}/output/{target2}_fit.root", "RECREATE")
    vtree = ROOT.TTree("vtree", "Tree holding fit parameters")
    # Create branches for fit parameters
    t_branch = np.zeros(1, dtype='float32')
    sigma_branch = np.zeros(1, dtype='float32')
    curve_branch = np.zeros(len(z), dtype='float32')  # Assumed length of curve

    vtree.Branch("t", t_branch, "t/F")
    vtree.Branch("sigma", sigma_branch, "sigma/F")
    vtree.Branch("curve", curve_branch, f"curve[{len(z)}]/F")  # Change length if needed

    t_branch[0] = bestfit_params_conv.t
    sigma_branch[0] = bestfit_params_conv.sigma    
    curve_branch[:] = bestfit_params_conv.curve
    vtree.Fill()

    # Fill the tree with data
    for param in params_list_conv:
        t_branch[0] = param.t        
        sigma_branch[0] = param.sigma        
        curve_branch[:] = param.curve
        vtree.Fill()
else:
    fit_file1 = ROOT.TFile(f"{beamtime}/bp-p-{target1}/output/{target1}_fit.root", "RECREATE")
    vtree = ROOT.TTree("vtree", "Tree holding fit parameters")
    # Create branches for fit parameters
    Phi0_branch = np.zeros(1, dtype='float32')
    R0_branch = np.zeros(1, dtype='float32')
    sigma_branch = np.zeros(1, dtype='float32')
    epsilon_branch = np.zeros(1, dtype='float32')
    curve_branch = np.zeros(len(z), dtype='float32')  # Assumed length of curve

    vtree.Branch("Phi0", Phi0_branch, "Phi0/F")
    vtree.Branch("R0", R0_branch, "R0/F")
    vtree.Branch("sigma", sigma_branch, "sigma/F")
    vtree.Branch("epsilon", epsilon_branch, "epsilon/F")
    vtree.Branch("curve", curve_branch, f"curve[{len(z)}]/F")  # Change length if needed

    Phi0_branch[0] = bestfit_params.Phi0
    R0_branch[0] = bestfit_params.R0
    sigma_branch[0] = bestfit_params.sigma
    epsilon_branch[0] = bestfit_params.epsilon
    curve_branch[:] = bestfit_params.curve
    vtree.Fill()

    # Fill the tree with data
    for param in params_list:
        Phi0_branch[0] = param.Phi0
        R0_branch[0] = param.R0
        sigma_branch[0] = param.sigma
        epsilon_branch[0] = param.epsilon
        curve_branch[:] = param.curve
        vtree.Fill()

plt.show()
if fit_file1.IsOpen():
    print("Fit file is successfully written.")
else:
    print("Failed to write the fit file.")
fit_file1.Write()
fit_file1.Close()