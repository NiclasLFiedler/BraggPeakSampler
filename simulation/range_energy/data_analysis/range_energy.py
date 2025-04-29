import ROOT
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uproot  # uproot is a great library for reading ROOT files in Python
matplotlib.use('TkAgg')  # or 'Qt5Agg'

class RangeEnergyRelationship:
    def __init__(self):
        # Initialize empty lists for energy, range, and standard deviation
        self.energy = []
        self.range = []
        self.std_dev = []

    def add_data(self, energy, range_value, std_dev):
        # Append the values to the corresponding lists
        self.energy.append(energy)
        self.range.append(range_value)
        self.std_dev.append(std_dev)
        
    def __repr__(self):
        # Provide a string representation of the stored data for easy visualization
        return "\n".join(rf"\SI{e:.3f}{{\mega\electronvolt}} & \SI{r:.3f}{{\milli\meter}} $\pm$ \SI{s:.3f}{{\milli\meter}}" 
                         for e, r, s in zip(self.energy, self.range, self.std_dev))

def range_energy_relationship(energy, alpha, p):
    return alpha * energy**p

# Define a Gaussian function
def gaussian(x, mean, sigma, amplitude):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

def get_range_energy(energy, path, enable_output=False, enable_plot=False):
    # File path to your ROOT file
    root_file_path = f'{path}/range_{energy}.root'
    
    # Open the ROOT file using uproot
    with uproot.open(root_file_path) as file:
        tree = file["range"]  # Replace 'range' with your tree name
    
        # Extract data from the 'range' branch
        values = tree["range"].array(library="np")  # Replace 'range' with your branch name
    
    # Calculate min and max values
    min_value = np.min(values)
    max_value = np.max(values)
    
    # Calculate range and determine the number of bins
    data_range = max_value - min_value
    if data_range < 4:
        data_range *= 20
    if data_range < 50:
        data_range *= 10
    n_bins = int(data_range * 10)  # Example: 10 bins per unit range
    
    # Create histogram data manually
    bin_contents, bin_edges = np.histogram(values, bins=n_bins, range=(min_value, max_value))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Find the peak value
    max_bin = np.argmax(bin_contents)
    peak_value = bin_centers[max_bin]

    initial_sigma_guess = np.std(values)
    # Fit a Gaussian to the histogram data
    initial_guess = [peak_value, initial_sigma_guess, np.max(bin_contents)]
    popt, pcov = curve_fit(gaussian, bin_centers, bin_contents, p0=initial_guess, maxfev=1000000)
    
    # Extract fit parameters and their standard deviations
    mean, sigma, amplitude = popt
    mean_stddev, sigma_stddev, amplitude_stddev = np.sqrt(np.diag(pcov))
    
    if enable_output:
        print(f"-------------Energy: {energy}")
        print(f"Initial guess: {initial_guess}")
        print()
        print(f"Peak value of the histogram: {peak_value:.3e}")
        print(f"Gaussian fit parameters:")
        print(f"Mean: {mean:.3e} ± {mean_stddev:.3e}")
        print(f"Sigma: {sigma:.3e} ± {sigma_stddev:.3e}")
        print(f"Amplitude: {amplitude:.3e} ± {amplitude_stddev:.3e}")
        print("--------------------------------------")
    
    # Plot the histogram and the Gaussian fit
    if enable_plot:
        plt.rcParams.update({'font.size': 15})
        plt.figure(figsize=(12, 6))
        plt.hist(bin_centers, bins=n_bins, weights=bin_contents, alpha=0.6, label='Data')
        plt.plot(bin_centers, gaussian(bin_centers, *popt), 'r-', label=f'Gaussian Fit\n' rf'$\mu$: {mean:.3e} mm ± {mean_stddev:.3e} mm' f'\n' rf'$\sigma$: {sigma:.3e} mm ± {sigma_stddev:.3e} mm')
        plt.xlabel('Range / mm')
        plt.ylabel('Counts')
        plt.title(f'{energy} MeV Proton Range with Gaussian Fit')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'{path}/range_fit_{energy}.pdf', format="pdf", bbox_inches="tight")
        plt.close()
        #plt.show()
    
    return (energy, mean, abs(sigma))

def analyse_range_energy(data : RangeEnergyRelationship, comparison_data: RangeEnergyRelationship = None):
   # Extract energies and ranges from the EnergyRangeData object
    energies = np.array(data.energy)
    ranges = np.array(data.range)
    
    # Perform curve fitting
    popt, pcov = curve_fit(range_energy_relationship, energies, ranges, maxfev = 1000000)
    
    # popt contains the optimized values for alpha and p
    alpha, p = popt
    std_devs = np.sqrt(np.diag(pcov))
    alpha_stddev, p_stddev = std_devs
    
    # Print the fitting parameters with their standard deviations
    print(f"Fitted parameters:")
    print(f"alpha = {alpha:.3e} ± {alpha_stddev:.3e}")
    print(f"p = {p:.3e} ± {p_stddev:.3e}")
    
    # Generate fitted ranges using the optimized parameters
    fitted_ranges = range_energy_relationship(energies, *popt)
    
    # Plot the original data
    plt.rcParams.update({'font.size': 20})
    plt.figure(figsize=(16, 8))
    plt.errorbar(energies, ranges, yerr=data.std_dev, fmt='o', ecolor='r', capsize=5, label=f'{name} Simulation')

    # Plot the fitted curve
    plt.plot(energies, fitted_ranges, 'b-', label=f'Fit {name}: $\\alpha_{{{name}}} = {alpha:.3e}$ $\\frac{{mm}}{{MeV^p}}$; $p_{{{name}}}$ = {p:.3e}')
 
    # If comparison data is provided, plot it
    if comparison_data:
        comp_energies = np.array(comparison_data.energy)
        comp_ranges = np.array(comparison_data.range)
        
        # Perform curve fitting for the comparison data
        comp_popt, comp_pcov = curve_fit(range_energy_relationship, comp_energies, comp_ranges, maxfev=1000000)
        
        comp_alpha, comp_p = comp_popt
        comp_std_devs = np.sqrt(np.diag(comp_pcov))
        comp_alpha_stddev, comp_p_stddev = comp_std_devs
        
        # Print the fitting parameters with their standard deviations
        print(f"comp_Fitted parameters:")
        print(f"alpha = {comp_alpha:.3e} ± {comp_alpha_stddev:.3e}")
        print(f"p = {comp_p:.3e} ± {comp_p_stddev:.3e}")

        # Generate fitted ranges for comparison data
        comp_fitted_ranges = range_energy_relationship(comp_energies, *comp_popt)
        
        # Plot the comparison data
        plt.errorbar(comp_energies, comp_ranges, yerr=comparison_data.std_dev, fmt='o', ecolor='g', capsize=5, label='ICRU Data')
        plt.plot(comp_energies, comp_fitted_ranges, label=f'Fit ICRU 49: $\\alpha_{{ICRU}} = {comp_popt[0]:.3e}$ $\\frac{{mm}}{{MeV^p}}$; $p_{{ICRU}}$ = {comp_popt[1]:.3e}')

    # Add labels, title, and legend
    plt.xlabel('Energy / MeV')
    plt.ylabel('Range / mm')
    plt.title(f'Range Energy: Simulated {name} and measured CSDA-Range H2O (ICRU 49)')
    plt.legend()

    # Show the plot
    plt.grid(True)
    plt.savefig(f"{path}/rangeenergy.pdf", format="pdf", bbox_inches="tight")
    #plt.savefig(f"air/rangeenergy.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
path = "h2o" #h2o pbwo4 air DSB EJ256
name = "H2O" #H2O PbWO4 AIR DSB EJ-256
icru_el = 0 # 0=h20, 1=air 

def main():
    data = RangeEnergyRelationship()
    energies = [3, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250]
    
    
    rho = []
    rho_h2o = 0.997
    rho_h2o_icru = 1.000
    rho_h2o_air=1.2048e-3
    rho.append(rho_h2o_icru)
    rho.append(rho_h2o_air)
    
    ICRU_ranges = []
    CSDA_ranges_h20 = [1.499e-2, 3.623e-2, 1.23e-1, 2.539e-1, 4.26e-1, 8.853e-1, 1.489, 2.227, 3.093, 4.080, 5.184, 6.398, 7.718, 1.146e1, 1.577e1, 2.062e1, 2.596e1, 3.174e1, 3.794e1] #g/cm^2
    CSDA_ranges_air= [1.737e-2, 4.173e-2, 1.408e-1, 2.899e-1, 4.855e-1, 1.007, 1.691, 2.528, 3.509, 4.628, 5.876, 7.250, 8.744, 1.297e1, 1.786e1, 2.334e1, 2.937e1, 3.590e1, 4.290e1]
    ICRU_ranges.append(CSDA_ranges_h20)
    ICRU_ranges.append(CSDA_ranges_air)
    
    icru_data= RangeEnergyRelationship()
    cutoff = 0
    # get_range_energy(15, enable_output=True, enable_plot=True)
    for index, energy in enumerate(energies):
        if(index == len(energies)-cutoff):
            break
        data.add_data(*get_range_energy(energy, path, enable_output=False, enable_plot=False))

    # for index, energy in enumerate(energies):
    #     print(rf"{energy} MeV ${data.range[index]:.3f} \pm {data.std_dev[index]:.3f}$")

    for index, energy in enumerate(energies):
        print(f"$\SI{{{energy}}}{{\mega\electronvolt}}$ & ${ICRU_ranges[icru_el][index]/rho[icru_el]*10:.2e} \pm {ICRU_ranges[icru_el][index]/rho[icru_el]*10*0.015:.2e}$ &")

    for index, CSDA_range in enumerate(ICRU_ranges[icru_el]):
        if(index == len(ICRU_ranges[icru_el])-cutoff):
            break
        icru_data.add_data(energies[index], CSDA_range/rho[icru_el]*10, CSDA_range/rho[icru_el]*0.015*10)

    analyse_range_energy(data, icru_data)

if __name__ == "__main__":
    main()