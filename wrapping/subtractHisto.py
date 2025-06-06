import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

def rebin_histogram(values, edges, rebin_factor):
    if len(values) % rebin_factor != 0:
        trim = len(values) % rebin_factor
        values = values[:-trim]
        edges = edges[:-(trim)]

    new_values = values.reshape(-1, rebin_factor).sum(axis=1)
    new_edges = edges[::rebin_factor]
    new_edges = np.append(new_edges, edges[-1])
    
    return new_values, new_edges

def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

def fit_uncertainty(x, A, tau, cov):
    dA2 = cov[0, 0]
    dtau2 = cov[1, 1]
    covAtau = cov[0, 1]

    df_dA = np.exp(-x / tau)
    df_dtau = A * x * np.exp(-x / tau) / tau**2

    sigma2 = (df_dA**2) * dA2 + (df_dtau**2) * dtau2 + 2 * df_dA * df_dtau * covAtau
    return np.sqrt(sigma2)

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def main():
    histname = 'h_ampl0'
    dir_path = '3005'
    measurement = "df2000"
    source = "sr90"
    output_path = f'{dir_path}/{measurement}_{source}.png'

    file1_path = f"{dir_path}/pwo_c0_{source}_{measurement}_300.root"
    file2_path = f"{dir_path}/pwo_c0_{measurement}_300.root"

    # Load histograms
    with uproot.open(file1_path) as file1, uproot.open(file2_path) as file2:
        h1 = file1[histname]
        h2 = file2[histname]

        values1, edges = h1.to_numpy()
        values2, _ = h2.to_numpy()

    values_diff = values1 - values2
    bin_centers = (edges[:-1] + edges[1:]) / 2
    
    if(source == "na22"):
        values_diff, edges = rebin_histogram(values_diff, bin_centers, 8)
        
        bin_centers = (edges[:-1] + edges[1:]) / 2

        
        A_guess = np.max(values_diff)
        mu_guess = bin_centers[np.argmax(values_diff)]
        sigma_guess = np.std(bin_centers)

        mask = values_diff > 0
        x_fit = bin_centers[mask]
        y_fit = values_diff[mask]

        try:
            popt, pcov = curve_fit(gauss, x_fit, y_fit, p0=[A_guess, mu_guess, sigma_guess])
            A_fit, mu_fit, sigma_fit = popt
            print(f"Fit successful: μ = {mu_fit:.2f}, σ = {sigma_fit:.2f}")
        except RuntimeError:
            print("Fit failed!")
            popt = None
    elif(source == "sr90"):
        nonzero_indices = np.nonzero(values_diff)[0]
        print("Nonzero indices:", len(nonzero_indices))
        
        total_counts = np.sum(values_diff)
        cum_counts = np.cumsum(values_diff)

        tail_threshold = 0.99  # Use last 15% of cumulative counts
        tail_start_index = np.searchsorted(cum_counts, tail_threshold * total_counts)
        
        x_tail = bin_centers[tail_start_index:]
        y_tail = values_diff[tail_start_index:]
        
        sigma_y = np.sqrt(np.maximum(y_tail, 1))        
        popt, pcov = curve_fit(exp_decay, x_tail, y_tail, p0=(1e7, 1000), bounds=([0, 10], [np.inf, 1e6]), sigma=sigma_y, 
        absolute_sigma=True, maxfev=1000000)        
        residuals = y_tail - exp_decay(x_tail, *popt)
        chi2 = np.sum((residuals / sigma_y)**2)
        ndf = len(x_tail) - len(popt)
        print(f"Chi²/ndf = {chi2:.2f}/{ndf} = {chi2/ndf:.2f}")
        perr = np.sqrt(np.diag(pcov))
        def exp_decay_upper(x, A, tau, dA, dtau):
            return (A + dA) * np.exp(-x / (tau - dtau))
        def exp_decay_lower(x, A, tau, dA, dtau):
            return (A - dA) * np.exp(-x / (tau + dtau))
        mask = exp_decay(x_tail, *popt) > 0.1
        y_upper = exp_decay_upper(x_tail[mask], *popt, perr[0], perr[1])
        y_lower = exp_decay_lower(x_tail[mask], *popt, perr[0], perr[1])
        # Calculate the intersection of the fit with y = 1
        A, tau = popt
        dA, dtau = perr
        # Intersection point (x value where y = 1)
        x_intersect = -tau * np.log(1 / A)
        # Propagate uncertainty
        dx_dA = tau / A  # Partial derivative of x with respect to A
        dx_dtau = -np.log(1 / A)  # Partial derivative of x with respect to tau
        x_intersect_std = np.sqrt((dx_dA * dA)**2 + (dx_dtau * dtau)**2)
        print(f"Intersection with y=1: x = {x_intersect:.2f} ± {x_intersect_std:.2f}")
        # Optionally, you can add this information to the plot
    
    # Plot
    plt.figure(figsize=(10,6))
    plt.step(edges[:-1], values_diff, where='mid', label='Subtracted Histogram', alpha=0.8)
    
    if source == "na22":
        x_dense = np.linspace(min(bin_centers), max(bin_centers), 1000)
        mask = gauss(x_dense, *popt) > 0.1
        plt.plot(x_dense[mask], gauss(x_dense, *popt)[mask], 'r-', label='Gaussian Fit')
    elif source == "sr90":
        x_dense = np.linspace(min(bin_centers), max(bin_centers), 1000)
        plt.plot(x_tail[mask], exp_decay(x_tail[mask], *popt), 'r--', label=f"Fit: A*exp(-x/τ)\nA={popt[0]:.2e}, τ={popt[1]:.2f}")
        plt.fill_between(x_tail[mask], y_lower, y_upper, color='red', alpha=0.3, label="Fit ± 1σ")
        plt.axvline(x_intersect, color='blue', linestyle='--', label=f"Intersection: x={x_intersect:.2f} ± {x_intersect_std:.2f}")


    plt.xlabel("X")
    plt.ylabel("Counts")
    plt.title("Histogram Subtraction and Gaussian Fit")
    plt.legend()
    plt.yscale('log')
    plt.grid(True)
    plt.tight_layout()
    plt.show()  # Keeps the window open until closed

if __name__ == "__main__":
    main()
