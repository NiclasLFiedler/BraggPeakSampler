import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

def range_energy_relationship(energy, alpha, p):
    return alpha * energy**p

def fit_range_peak(ranges, bins=200):

    hist, edges = np.histogram(ranges, bins=bins)

    centers = 0.5*(edges[:-1]+edges[1:])
    peak_idx = np.argmax(hist)

    peak_pos = centers[peak_idx]


    # estimate sigma from distribution width
    sigma_guess = np.std(ranges)

    # fit window around peak
    fit_range = 3*sigma_guess

    mask = (centers > peak_pos-fit_range) & (centers < peak_pos+fit_range)


    p0 = [hist[peak_idx], peak_pos, sigma_guess]


    popt, pcov = curve_fit(gauss, centers[mask], hist[mask], p0=p0)
    perr = np.sqrt(np.diag(pcov))

    A, mu, sigma = popt

    return mu, sigma, A, centers, hist, popt, perr

materials = ["pbwo4proj", "h2oproj"]
# materials = ["pbwo4", "h2o"]

energies = [20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200]

useMaterial = 1
material = materials[useMaterial]

files = [f"{material}/raw_data_{energy}.root" for energy in energies]

meanRanges = []

for file in files:
    tree = uproot.open(file)["braggsampler"]

    df = tree.arrays(
        ["event","pos","eDep", "dEdX", "trackid","eKin","eTot", "StepLength"],
        library="pd"
    )

    # only primary proton
    df = df[df.trackid == 1]
    
    ranges = []
    ranges_path = []

    for event, g in df.groupby("event"):

        g = g.sort_values("pos")

        zend = g.pos.iloc[-1]

        ranges.append(zend)
        ranges_path.append(g.StepLength.sum())

    ranges = np.array(ranges)
    ranges_path = np.array(ranges_path)

    fit_mu, fit_sigma, fit_A, centers, hist, popt, perr = fit_range_peak(ranges, 1000)
    fit_mu_path, fit_sigma_path, fit_A_path, centers_path, hist_path, popt_path, perr_path = fit_range_peak(ranges_path, 1000)

    meanRanges.append((fit_mu, fit_sigma, fit_mu_path, fit_sigma_path, perr))

    print("Mean projected range =", fit_mu)
    print("Std =", fit_sigma)
    print("Mean path length =", fit_mu_path)
    print("Std =", fit_sigma_path)

    # plt.figure(figsize=(12, 5))
    # plt.subplot(1, 2, 1)
    # plt.bar(centers, hist, width=centers[1]-centers[0], alpha=0.5, label="Histogram")
    # plt.plot(centers, gauss(centers, *popt), color='red', label='Gaussian Fit')
    # plt.title(f"Projected Range Distribution for {material} at {220} MeV")
    # plt.xlabel("Projected Range (mm)")
    # plt.ylabel("Counts")
    # plt.legend()
    # plt.show()

CSDAranges=  [x[0] for x in meanRanges]
popt, pcov = curve_fit(range_energy_relationship, energies, CSDAranges, maxfev = 1000000)

alpha, p = popt
std_devs = np.sqrt(np.diag(pcov))
alpha_stddev, p_stddev = std_devs

print(f"\nPbwo4 Fitted parameters:")
print(f"alpha = {alpha:.3e} ± {alpha_stddev:.3e}")
print(f"p = {p:.3e} ± {p_stddev:.3e}")     
print(f"Covariance Matrix:{pcov[0][1]}")


plt.plot(energies, CSDAranges, 'o', label="Projected Range")
plt.plot(energies, range_energy_relationship(energies, alpha, p), label="Fitted Curve")
plt.xlabel("Energy (MeV)")
plt.ylabel("Projected Range (mm)")
plt.legend()
plt.show()
