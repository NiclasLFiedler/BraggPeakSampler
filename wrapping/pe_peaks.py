import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def gaussian(x, a, mu, sigma):
    return a * np.exp(-(x - mu)**2 / (2 * sigma**2))

def load_histogram(root_file, hist_name):
    with uproot.open(root_file) as file:
        hist = file[hist_name]
        counts, edges = hist.to_numpy()
        centers = 0.5 * (edges[:-1] + edges[1:])
        return centers, counts

def fit_peaks(x, y, peak_indices, window=10, nbPeaks=5):
    fits = []
    for idx in peak_indices[:nbPeaks]:  # Limit to first nbPeaks
        center = x[idx]
    
        mask = (x > center - window) & (x < center + window)
        x_fit = x[mask]
        y_fit = y[mask]

        try:
            p0 = [y[idx], center, 1.0]
            popt, _ = curve_fit(gaussian, x_fit, y_fit, p0=p0)
            fits.append(popt)
        except RuntimeError:
            print(f"Fit failed for peak at index {idx}")
            continue
    return fits

def plot_results(x, y, fits):
    plt.figure(figsize=(10, 6))
    plt.step(x, y, label="Data: PE Linear Fit:")
    plt.plot(x[peak_indices], y[peak_indices], 'rx', label="Found   Peaks", markersize=8, markeredgewidth=2)
    
    
    x_dense = np.linspace(x[0], x[-1], 100000)
    for i, (a, mu, sigma) in enumerate(fits):
        mask = gaussian(x_dense, a, mu, sigma) > 0.1
        plt.plot(x_dense[mask], gaussian(x_dense, a, mu, sigma)[mask], color='orange', label=rf"Peak {i+1}: $\mu$ = {mu:.2f}; $\sigma$ = {sigma:.2f}")
        plt.axvline(mu, color='gray', linestyle='--', alpha=0.5)

    plt.xlabel("ADC / Channels")
    plt.ylabel("Counts")
    plt.title("SiPM Signal: Photoelectron Peaks")
    plt.legend()
    plt.yscale('log')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_linear_fit():
    plt.figure(figsize=(8, 5))
    plt.errorbar(pe_numbers, mu_positions, yerr=sigma_mu, fmt='o', capsize=5, color='blue', label=rf'$\mu \pm \sigma$ from Gauss fits')
    plt.plot(pe_numbers, slope * pe_numbers + intercept, 'r-', label=f'Linear fit:\n' fr'$\mu$ = {slope:.2f}·PE {intercept:+.2f}')
    plt.xlabel('Photoelectron Number')
    plt.ylabel('μ Position (ADC Channels)')
    plt.title('Photoelectron Calibration')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    nbPeaks = 6
    peakWidthEst = 5 
    peakDistanceEst = 5

    hist_dir = '3005' 
    root_file = f"{hist_dir}/no_crystal_20_led_120s.root"
    hist_name = f"h_ampl0"

    x, y = load_histogram(root_file, hist_name)


    peak_indices, _ = find_peaks(y, width=peakWidthEst, distance=peakDistanceEst)

    fits = fit_peaks(x, y, peak_indices, peakWidthEst, nbPeaks)

    mu_positions = [mu for _, mu, _ in fits]
    sigma_mu = np.array([sigma for _, _, sigma in fits])
    pe_numbers = np.arange(1, len(mu_positions) + 1)
    slope, intercept = np.polyfit(pe_numbers, mu_positions, 1)

    plot_results(x, y, fits)
    plot_linear_fit()
