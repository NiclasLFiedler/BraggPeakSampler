import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
from scipy.stats import exponnorm
from scipy.stats import moyal
from dataclasses import dataclass
import sys
sys.path.append("../../range_energy/data_analysis")
import analysisFunctions

@dataclass
class TargetParameters:
    Thickness: float
    Pmod_theo: float
    range: float
    sigma: float
    sigma_t: float
    t: float
    resRange: float
    sigma_res: float
    sigma_t_res: float
    t_res: float
    Pmod_res: float
    Pmod_sim: float
    energy: float
    sigma_E: float
    sigma_T_E: float
    t_E: float
    Pmod_E: float

def print_target_parameters(targetParamsList):

    print("\n" + "=" * 190)
    print("Target Parameters")
    print("=" * 190)

    print(
        f"{'Thickness':>10} "
        f"{'Pmod_theo':>10} "
        f"{'Range':>10} "
        f"{'ResRange':>10} "
        f"{'Sigma':>10} "
        f"{'Sigma_res':>10} "
        f"{'t':>10} "
        f"{'t_res':>10} "
        f"{'Sigma_t':>10} "
        f"{'Sigma_t_res':>10} "
        f"{'Pmod_res':>10} "
        f"{'Pmod_sim':>10} "
        f"{'Energy':>10} "
        f"{'Sigma_E':>10} "
        f"{'Sigma_T_E':>10} "
        f"{'t_E':>10} "
        f"{'Pmod_E':>10}"
    )

    print("-" * 190)

    for target in targetParamsList:
        print(
            f"{target.Thickness:10.1f} "
            f"{target.Pmod_theo:10.1f} "
            f"{target.range:10.3f} "
            f"{target.resRange:10.3f} "
            f"{target.sigma:10.3f} "
            f"{target.sigma_res:10.3f} "
            f"{target.t:10.3f} "
            f"{target.t_res:10.3f} "
            f"{target.sigma_t_res:10.3f} "
            f"{target.sigma_t:10.3f} "
            f"{target.Pmod_res:10.3f} "
            f"{target.Pmod_sim:10.3f} "
            f"{target.energy:10.3f} "
            f"{target.sigma_E:10.3f} "
            f"{target.sigma_T_E:10.3f} "
            f"{target.t_E:10.3f} "
            f"{target.Pmod_E:10.3f}"
        )

    print("=" * 190)
    
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
    E = np.asarray(E_MeV)
    range_cm = np.zeros_like(E, dtype=float)
    mask_valid = E >= 1.0

    data = np.load("../../range_energy/data_analysis/h2o_range_energy.npz")
    useSumFit = data["useSumFit"]

    if useSumFit:
        alpha = data["alpha"]
        if np.any(mask_valid):
            E_valid = E[mask_valid] if isinstance(E, np.ndarray) else E
            range_cm = analysisFunctions.range_energy_sum(E_valid, *alpha)
    else:
        alpha = data["alpha"]
        p = data["p"]
        if np.any(mask_valid):
            E_valid = E[mask_valid] if isinstance(E, np.ndarray) else E
            range_cm = analysisFunctions.range_energy_relationship(E_valid, alpha[0], p)
        
    return range_cm
 
 
def gaussian(x, A, mu, sigma):
    """Gaussian function for fitting."""
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
 
 
def _fit_gaussian_to_data(data, bounds_upper=None):
    """
    Fit a Gaussian to histogram data.
    
    Args:
        data: numpy array of data points
        bounds_upper: optional upper bounds for [A, mu, sigma]
    
    Returns:
        popt: fitted parameters [A, mu, sigma]
        or (np.nan, np.nan, np.nan) if fit fails
    """
    try:
        mu = np.mean(data)
        sigma = np.std(data)
        A = 1 / (sigma * np.sqrt(2 * np.pi))
        
        hist, bin_edges = np.histogram(data, bins=4000, density=True)
        bincenters = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        
        if bounds_upper is None:
            bounds_upper = [np.inf, np.inf, np.inf]
        
        popt, pcov = curve_fit(
            gaussian, 
            bincenters, 
            hist,
            p0=[A, mu, sigma],
            bounds=([0, 0, 0], bounds_upper),
            maxfev=5000
        )
        return popt
    except Exception as e:
        print(f"Warning: Could not fit Gaussian: {e}")
        return (np.nan, np.nan, np.nan)


def _plot_histogram_with_fit(ax, data, popt, xlabel, ylabel, title, color, unit):
    """
    Plot histogram with Gaussian fit overlay.
    
    Args:
        ax: matplotlib axis
        data: numpy array of data points
        popt: fitted parameters [A, mu, sigma]
        xlabel: label for x-axis
        ylabel: label for y-axis
        title: plot title
        color: histogram color
        unit: unit string for the data (e.g., 'MeV', 'cm', 'mm')
    """
    # Create histogram
    ax.hist(
        data, 
        bins=4000, 
        density=True, 
        alpha=0.7, 
        color=color, 
        edgecolor='black',
        label='Data',
        histtype='step'
    )
    
    # Plot fit if valid
    if not np.isnan(popt[0]):
        x = np.linspace(data.min(), data.max(), 4000)
        y = gaussian(x, *popt)
        ax.plot(x, y, 'r-', linewidth=2, label='Gaussian Fit')
        
        # Add text box with fit parameters
        ax.text(0.95, 0.95, 
                f'μ = {popt[1]:.3f} {unit}\nσ = {popt[2]/popt[1]*100:.3f} %',
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)


def _print_summary_statistics(data_dict):
    """
    Print summary statistics for all datasets.
    
    Args:
        data_dict: dict with keys like 'Kinetic Energy (MeV)', values are tuples of (data, popt)
    """
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    for name, (data, popt, unit) in data_dict.items():
        print(f"\n{name}:")
        print(f"  Mean:   {np.mean(data):.4f}")
        print(f"  Median: {np.median(data):.4f}")
        print(f"  Std:    {np.std(data):.4f}")
        print(f"  Min:    {np.min(data):.4f}")
        print(f"  Max:    {np.max(data):.4f}")
        
        if not np.isnan(popt[0]):
            print(f"\n{name} Gaussian Fit:")
            print(f"  μ = {popt[1]:.4f} {unit}")
            print(f"  σ = {popt[2]:.4f} {unit}")
            print(f"  FWHM = {2.355 * popt[2]:.4f} {unit}")
    
    print("="*60)


def analyse_data(file_path, enablePrint=False, enablePlot=False):
    """
    Analyze GEANT4 proton data from ROOT file.
    
    Args:
        file_path: path to ROOT file
        enablePrint: whether to print summary statistics
        enablePlot: whether to plot histograms
    
    Returns:
        fitParameters: list of fitted parameters for [eKin, range, resRange]
    """
    fitParameters = []
    
    # Read data from ROOT file
    try:
        file = uproot.open(file_path)
        tree = file["braggsampler"]
        eKin = tree["eKin"].array(library="np")
        resRange = tree["CSDARange"].array(library="np")
        
        if enablePrint:
            print(f"Successfully read {len(eKin)} events")
            print(f"Kinetic energy range: {np.min(eKin):.3f} - {np.max(eKin):.3f} MeV")
            print(f"Residual range range: {np.min(resRange):.3f} - {np.max(resRange):.3f} mm")
        
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
    
    # Fit Gaussians to all datasets
    popt_eKin = _fit_gaussian_to_data(eKin, bounds_upper=[np.inf, 500, np.inf])
    popt_range = _fit_gaussian_to_data(range_cm, bounds_upper=[np.inf, 500, np.inf])
    popt_resRange = _fit_gaussian_to_data(resRange, bounds_upper=[np.inf, 500, np.inf])
    
    fitParameters = [popt_eKin, popt_range, popt_resRange]
    
    # Create and populate plots if enabled
    if enablePlot:
        fig, axes = plt.subplots(1, 3, figsize=(21, 5))
        
        _plot_histogram_with_fit(
            axes[0], eKin, popt_eKin,
            xlabel='Kinetic Energy (MeV)',
            ylabel='Probability Density',
            title='Kinetic Energy Distribution',
            color='blue',
            unit='MeV'
        )
        
        _plot_histogram_with_fit(
            axes[1], range_cm, popt_range,
            xlabel='Range in Water (cm)',
            ylabel='Probability Density',
            title='Range Distribution (Bortfeld)',
            color='green',
            unit='cm'
        )
        
        _plot_histogram_with_fit(
            axes[2], resRange, popt_resRange,
            xlabel='Residual Range (mm)',
            ylabel='Probability Density',
            title='Residual Range Distribution',
            color='purple',
            unit='mm'
        )
        
        plt.tight_layout()
        plt.show()
    else:
        plt.close()
    
    # Print summary statistics if enabled
    if enablePrint:
        data_dict = {
            'Kinetic Energy (MeV)': (eKin, popt_eKin, 'MeV'),
            'Range in Water (cm)': (range_cm, popt_range, 'cm'),
            'Residual Range (mm)': (resRange, popt_resRange, 'mm'),
        }
        _print_summary_statistics(data_dict)
    
    return fitParameters


def main():
    data_file_path = "data"
    data_file_ref = "data_0_0.root"

    filepathref = f"{data_file_path}/{data_file_ref}"
 
    targetThicknesses = range(100, 150, 100)
    pmods = range(100, 900, 200)
 
    targetParamsList = []
 
    # Reference measurements
    fitparams = analyse_data(filepathref, True, True)
    targetParamsList.append(TargetParameters(
        Thickness=0,
        Pmod_theo=0,
        range=fitparams[1][1],
        sigma=fitparams[1][2],
        sigma_t=0,
        t=0,
        Pmod_sim=0,
        resRange=fitparams[2][1],
        sigma_res=fitparams[2][2],
        sigma_t_res=0,
        t_res=0,
        Pmod_res=0,
        energy=fitparams[0][1],
        sigma_E=fitparams[0][2],
        sigma_T_E=0,
        t_E=0,
        Pmod_E=0
    ))
 
    # Scan over thicknesses and modulation powers
    for thickness in targetThicknesses:
        for pmod in pmods:   
            dataFile = f"{data_file_path}/data_{pmod}_{thickness}.root"         
            fitparams = analyse_data(dataFile)
 
            # Range analysis
            t = targetParamsList[0].range - fitparams[1][1]
            variance = fitparams[1][2]**2 - targetParamsList[0].sigma**2
            modulationPower = variance / t * 10000
 
            # Energy analysis
            tE = targetParamsList[0].energy - fitparams[0][1]
            varianceE = fitparams[0][2]**2 - targetParamsList[0].sigma_E**2
            EnergyModulation = varianceE / tE * 1000
 
            # Residual Range analysis
            t_res = targetParamsList[0].resRange - fitparams[2][1]
            variance_res = fitparams[2][2]**2 - targetParamsList[0].sigma_res**2
            modulationPower_resRange = variance_res / t_res * 10000
 
            targetParamsList.append(TargetParameters(
                Thickness=thickness,
                Pmod_theo=pmod,
                range=fitparams[1][1],
                sigma=fitparams[1][2],
                sigma_t=np.sqrt(variance),
                t=t,
                Pmod_sim=modulationPower,
                resRange=fitparams[2][1],
                sigma_res=fitparams[2][2],
                sigma_t_res=np.sqrt(variance_res),
                t_res=t_res,
                Pmod_res=modulationPower_resRange,
                energy=fitparams[0][1],
                sigma_E=fitparams[0][2],
                sigma_T_E=np.sqrt(varianceE),
                t_E=tE,
                Pmod_E=EnergyModulation
            ))
 
    print_target_parameters(targetParamsList)

 
if __name__ == "__main__":
    main()
 

