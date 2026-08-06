import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
from scipy.stats import exponnorm
from scipy.stats import moyal


def landau(x, A, mpv, sigma):
    """
    Landau distribution approximation using the Moyal distribution.

    Parameters
    ----------
    x : array
        x values
    A : float
        amplitude
    mpv : float
        most probable value
    sigma : float
        width parameter

    Returns
    -------
    y : array
        Landau-like distribution
    """
    return A * moyal.pdf(x, loc=mpv, scale=sigma)
 
def shiftedLandau(x, x0, A, mpv, sigma):
    return x0-landau(x, A, mpv, sigma)

def bortfeld_range(E_MeV):
    """
    Calculate range in water for proton kinetic energy using Bortfeld's formula.
    
    Parameters:
    E_MeV : float or array
        Kinetic energy in MeV
    
    Returns:
    range : float or array
        Range in cm (in water)
    
    Note: This is valid for E >= 1 MeV
    """
    E = np.asarray(E_MeV)
    
    # Bortfeld parameterization for protons in water
    # Valid range: roughly 1 MeV to 250 MeV
    range_cm = np.zeros_like(E, dtype=float)
    
    mask_valid = E >= 1.0
    
    # For E >= 1 MeV
    if np.any(mask_valid):
        E_valid = E[mask_valid] if isinstance(E, np.ndarray) else E

        p = 1.77
        alpha = 0.0225
        range_cm = alpha * np.power(E_valid, p)
        
        
        # if E_valid <= 2.5:
        #     R = 0.56 * E_valid - 0.01
        # else:
        #     R = 0.31 * E_valid**1.5 + 0.06
        
        # if isinstance(E, np.ndarray):
        #     range_cm[mask_valid] = R
        # else:
        #     range_cm = R

    
    return range_cm
 
 
def gaussian(x, A, mu, sigma):
    """Gaussian function for fitting."""
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
 
 
# ============================================================================
# Main Analysis
# ============================================================================
 
def analyse_data(file_path, enablePrint = False, enablePlot = False):

    fitParameters = []

    try:
        file = uproot.open(file_path)
        tree = file["braggsampler"]
        
        # Extract data
        eKin = tree["eKin"].array(library="np")
        
        if enablePrint:
            print(f"Successfully read {len(eKin)} events")
            print(f"Kinetic energy range: {np.min(eKin):.3f} - {np.max(eKin):.3f} MeV")
        
    except Exception as e:
        print(f"Error reading ROOT file: {e}")
        print("Make sure the file path is correct and the ROOT file is properly formatted.")
        return
    
    # Calculate range from kinetic energy using Bortfeld's formula
    if enablePrint:
        print("\nCalculating range from kinetic energy using Bortfeld's formula...")
    range_cm = bortfeld_range(eKin)
    
    if enablePrint:
        print(f"Range values: {np.min(range_cm):.3f} - {np.max(range_cm):.3f} cm")
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # ========================================================================
    # HISTOGRAM 1: Kinetic Energy with Gaussian Fit
    # ========================================================================
    ax1 = axes[0]
    
    # Create histogram
    counts_eKin, bins_eKin, patches = ax1.hist(
        eKin, 
        bins=2000, 
        density=True, 
        alpha=0.7, 
        color='blue', 
        edgecolor='black',
        label='Data',
        histtype='step'
    )
    
    # Fit Gaussian to kinetic energy
    try:
        # Estimate initial parameters
        mu_eKin = np.mean(eKin)
        sigma_eKin = np.std(eKin)
        A_eKin = 1 / (sigma_eKin * np.sqrt(2 * np.pi))

        hist, bin_edges = np.histogram(eKin, bins=2000, density=True)
        bincenters = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # Fit
        popt_eKin, pcov_eKin = curve_fit(
            gaussian, 
            bincenters, 
            hist,
            p0=[A_eKin, mu_eKin, sigma_eKin],
            maxfev=5000
        )
        fitParameters.append(popt_eKin)
        if enablePlot:
            x_eKin = np.linspace(eKin.min(), eKin.max(), 2000)
            y_eKin = gaussian(x_eKin, *popt_eKin)

            ax1.plot(x_eKin, y_eKin, 'r-', linewidth=2, label='Gaussian Fit')

            # Display fit parameters
            chi2_eKin = np.sum((np.histogram(eKin, bins=2000, density=True)[0] - 
                                gaussian(np.histogram(eKin, bins=2000)[1][:-1], *popt_eKin)) ** 2)

            ax1.text(0.95, 0.95, 
                    f'μ = {popt_eKin[1]:.3f} MeV\nσ = {popt_eKin[2]/popt_eKin[1]*100:.3f} %',
                    transform=ax1.transAxes,
                    fontsize=11,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
    except Exception as e:
        fitParameters.append((np.nan, np.nan, np.nan))
        print(f"Warning: Could not fit Gaussian to kinetic energy: {e}")

    if enablePlot:
        ax1.set_xlabel('Kinetic Energy (MeV)', fontsize=12)
        ax1.set_ylabel('Probability Density', fontsize=12)
        ax1.set_title('Kinetic Energy Distribution', fontsize=13, fontweight='bold')
        ax1.legend(fontsize=11)
        ax1.grid(True, alpha=0.3)
    
    # ========================================================================
    # HISTOGRAM 2: Range with Gaussian Fit
    # ========================================================================
    ax2 = axes[1]
    
    # Create histogram
    counts_range, bins_range, patches = ax2.hist(
        range_cm, 
        bins=2000, 
        density=True, 
        alpha=0.7, 
        color='green', 
        edgecolor='black',
        label='Data',
        histtype='step'
    )
    
    # Fit Gaussian to range
    try:
        # Estimate initial parameters
        mu_range = np.mean(range_cm)
        sigma_range = np.std(range_cm)
        A_range = 1 / (sigma_range * np.sqrt(2 * np.pi))
        
        hist, bin_edges = np.histogram(range_cm, bins=2000, density=True)
        bincenters = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        popt_range, pcov_range = curve_fit(
            gaussian, 
            bincenters, 
            hist,
            p0=[A_range, mu_range, sigma_range],
            maxfev=5000
        )
        fitParameters.append(popt_range)
        if enablePlot:
            x_range = np.linspace(range_cm.min(), range_cm.max(), 2000)
            y_range = gaussian(x_range, *popt_range)
            ax2.plot(x_range, y_range, 'r-', linewidth=2, label='Gaussian Fit')

            # Display fit parameters
            chi2_range = np.sum((np.histogram(range_cm, bins=2000, density=True)[0] - 
                                 gaussian(np.histogram(range_cm, bins=2000)[1][:-1], *popt_range)) ** 2)

            ax2.text(0.95, 0.95, 
                    f'μ = {popt_range[1]:.3f} cm\nσ = {popt_range[2]/popt_range[1]*100:.3f} %',
                    transform=ax2.transAxes,
                    fontsize=11,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
    except Exception as e:
        fitParameters.append((np.nan, np.nan, np.nan))
        print(f"Warning: Could not fit Gaussian to range: {e}")
    if enablePlot:
        ax2.set_xlabel('Range in Water (cm)', fontsize=12)
        ax2.set_ylabel('Probability Density', fontsize=12)
        ax2.set_title('Range Distribution (Bortfeld)', fontsize=13, fontweight='bold')
        ax2.legend(fontsize=11)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()
    else:
        plt.close()
    
    # Print summary statistics
    if enablePrint:
        print("\n" + "="*60)
        print("SUMMARY STATISTICS")
        print("="*60)
        print(f"\nKinetic Energy (MeV):")
        print(f"  Mean:   {np.mean(eKin):.4f}")
        print(f"  Median: {np.median(eKin):.4f}")
        print(f"  Std:    {np.std(eKin):.4f}")
        print(f"  Min:    {np.min(eKin):.4f}")
        print(f"  Max:    {np.max(eKin):.4f}")

        print(f"\nRange in Water (cm):")
        print(f"  Mean:   {np.mean(range_cm):.4f}")
        print(f"  Median: {np.median(range_cm):.4f}")
        print(f"  Std:    {np.std(range_cm):.4f}")
        print(f"  Min:    {np.min(range_cm):.4f}")
        print(f"  Max:    {np.max(range_cm):.4f}")

        if 'popt_eKin' in locals():
            print(f"\nKinetic Energy Gaussian Fit:")
            print(f"  μ = {popt_eKin[1]:.4f} MeV")
            print(f"  σ = {popt_eKin[2]:.4f} MeV")
            print(f"  FWHM = {2.355 * popt_eKin[2]:.4f} MeV")

        if 'popt_range' in locals():
            print(f"\nRange Gaussian Fit:")
            print(f"  μ = {popt_range[1]:.4f} cm")
            print(f"  σ = {popt_range[2]:.4f} cm")
            print(f"  FWHM = {2.355 * popt_range[2]:.4f} cm")

        print("="*60)

    return fitParameters

def main():
    data_file_path = "data"
    data_file_name = "data_0.root"
    data_file_name2 = "data_200.root"

    filepath = f"{data_file_path}/{data_file_name}"

    fitparams = []
    fitparams.append(analyse_data(filepath))
    filepath = f"{data_file_path}/{data_file_name2}"
    fitparams.append(analyse_data(filepath))

    print(fitparams[0][1])
    print()
    print(print(fitparams[1][1]))
    t = fitparams[0][1][1]-fitparams[1][1][1]

    print(t)
    variance = fitparams[1][1][2]**2-fitparams[0][1][2]**2
    print(variance)
    modulationPower = variance/t
    print(f"Modulation Power: {modulationPower*1000:.2f} um")
 
 
if __name__ == "__main__":
    main()
 

