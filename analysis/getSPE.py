import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import uproot

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

show = False
file = uproot.open("../data/paperBeamtime/rawData/SPEmeasurment.root")
tree = file["RawData"]
charge = tree["Charge"].array(library="np")
channel = tree["Channel"].array(library="np")
board = tree["Board"].array(library="np")   

eff_channel = channel + 16 * board

unique_channels = np.unique(eff_channel)

charges_per_channel = {ch: charge[eff_channel == ch] for ch in unique_channels}

histograms = {}      # store hist arrays
bin_edges = {}        # store bin edges

for ch in unique_channels:
    hist, bins = np.histogram(charges_per_channel[ch], bins=8000, range=(0, 32000))
    histograms[ch] = hist
    bin_edges[ch] = bins

# for ch in unique_channels:
#     bin_centers = 0.5 * (bin_edges[ch][1:] + bin_edges[ch][:-1])
#     entries = len(charges_per_channel[ch])

#     plt.step(bin_centers, histograms[ch], where="mid",label=f"Channel {ch} — {entries} entries")
#     plt.title(f"Channel {ch} Dark Spectrum")
#     plt.xlabel("Charge")
#     plt.yscale("log")
#     plt.grid(True)
#     plt.ylabel("Counts")
#     plt.legend()
#     plt.show()

# exit()

params = []

for selectCH in unique_channels:
    if selectCH == 8 or selectCH == 9 or selectCH == 0:
        continue

    bin_centers = 0.5 * (bin_edges[selectCH][1:] + bin_edges[selectCH][:-1])
    y = histograms[selectCH]

    distance=100
    width = 5

    peaks, properties = find_peaks(y, distance=distance, width = width)
    # plt.plot(bin_centers[peaks], y[peaks], "x", color="red", label="Detected Peaks")
    # plt.step(bin_centers, y, 'k', label='Histogram')
    # plt.grid(True)
    # plt.yscale("log")
    # plt.title(f"SiPM dark spectrum Channel {selectCH}")
    # plt.show()
    mask = (peaks > 50)
    peaks = peaks[mask]
    first_three = np.sort(peaks)[:3]
    fits = []

    for p in first_three:
        # Choose a fitting window around each peak
        left = max(0, p - 25)
        right = min(len(bin_centers), p + 25)

        x_fit = bin_centers[left:right]
        y_fit = y[left:right]

        # initial guesses
        A0 = y[p]
        mu0 = bin_centers[p]
        sigma0 = (bin_centers[1] - bin_centers[0]) * 5

        popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[A0, mu0, sigma0], maxfev=100000)

        fits.append((popt, np.sqrt(np.diag(pcov))))

    print(f"=== SPE Peak Fit Results {selectCH} ===")
    for i, ((A, mu, sigma), (dA, dmu, dsigma)) in enumerate(fits, start=1):
        print(f"Peak {i}:")
        print(f"   Mean μ{i}       = {mu:.3f} ± {dmu:.3f}")
        print(f"   Width σ{i}      = {sigma:.3f} ± {dsigma:.3f}")
        print("")

    # Extract μ and σ values
    mus = np.array([f[0][1] for f in fits])
    sigmas = np.array([f[0][2] for f in fits])

    mus_err = np.array([f[1][1] for f in fits])
    sigmas_err = np.array([f[1][2] for f in fits])

    # Gains between peaks

    pxt = 0.2
    print(f"Sigma 1: {sigmas[0]}")
    print(f"Sigma 2: {sigmas[1]}")
    sigma_gain = np.sqrt((sigmas[1]**2 - sigmas[0]**2))
    sigma_elec = np.sqrt(sigmas[0]**2-sigma_gain**2)

    print(f"Sigma gain {sigma_gain}")
    print(f"Sigma elec {sigma_elec}")
    if np.isnan(sigma_elec):
        print("ERROR: sigma_elec is NaN → stopping.")
        raise SystemExit
    gains = mus[1:] - mus[:-1]
    gains_err = np.sqrt(mus_err[1:]**2 + mus_err[:-1]**2)

    # Gain widths
    gain_sigmas = np.sqrt(sigmas[1:]**2 + sigmas[:-1]**2)
    gain_sigmas_err = np.sqrt(
        (sigmas[1:] * sigmas_err[1:] / gain_sigmas)**2 +
        (sigmas[:-1] * sigmas_err[:-1] / gain_sigmas)**2
    )

    params.append([])

    print("=== SiPM Physical Parameters ===")
    for i in range(len(gains)):
        print(f"Gain {i+1} (μ{i+2} – μ{i+1})       = {gains[i]:.3f} ± {gains_err[i]:.3f}")
        print(f"Gain width σ_gain{i+1}             = {gain_sigmas[i]:.3f} ± {gain_sigmas_err[i]:.3f}")
        print("")

    plt.figure(figsize=(9,5))
    plt.plot(bin_centers[peaks], y[peaks], "x", color="red", label="Detected Peaks")
    plt.step(bin_centers, y, 'k', label='Histogram')

    x_dense = np.linspace(bin_centers[0], bin_centers[-1], 2000)

    for popt, _ in fits:
        plt.plot(x_dense, gaussian(x_dense, *popt), label=f"Peak at μ={popt[1]:.2f}")

    plt.legend()
    plt.xlabel("ADC / Charge / Amplitude")
    plt.ylabel("Counts")
    plt.grid(True)
    plt.yscale("log")
    plt.ylim(0.1, max(y)*1.2)
    plt.title(f"SiPM dark spectrum with fitted SPE peaks Channel {selectCH}")
    if show:
        plt.show()
    else:
        plt.close()
