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
    bin_contents[0] = 0
    # Find the peak value
    max_bin = np.argmax(bin_contents[10:])
    peak_value = bin_centers[max_bin+10]

    initial_sigma_guess = np.std(values)
    # Fit a Gaussian to the histogram data
    initial_guess = [peak_value, initial_sigma_guess, np.max(bin_contents[10:])]
    print(initial_guess)
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
        plt.show()
    
    return (energy, mean, abs(sigma))

def GetWET(data, water_data):
    WET = []
    
    print(f"WaterRange {water_data[0]}")
    print()
    for idx in range(0,27):
        waterRange = water_data[0][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(27,28):
        waterRange = water_data[1][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(28,29):
        waterRange = water_data[2][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(29,30):
        waterRange = water_data[3][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(30,31):
        waterRange = water_data[4][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(31,32):
        waterRange = water_data[5][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
        
    WETDepth = []
    for idx, value in enumerate(WET):
        if (idx == 0):
            WETDepth.append(WET[idx]/2)
        else:
            WETDepth.append((WET[idx-1]+WET[idx])/2) 

    print("Water equivalent thicknesses")
    for idx, value in enumerate(WET):
        if idx > 0:
            print(f"CH: {idx+1}; WET: {WET[idx]-WET[idx-1]}")
        else:
            print(f"CH: {idx+1}; WET: {WET[idx]}")
    
    print(f"{WET}")
    print("Water equivalent Depth")
    print(WETDepth)


def GetWETtarget(data, water_data):
    WET = []
    
    print(f"WaterRange {water_data[0]}")
    print()
    for idx in range(0,27):
        waterRange = water_data[0][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(21,22):
        waterRange = water_data[1][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(22,23):
        waterRange = water_data[2][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(23,24):
        waterRange = water_data[3][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(24,25):
        waterRange = water_data[4][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
    print()
    for idx in range(25,26):
        waterRange = water_data[5][1]
        print(f"{data[idx][1]}")
        WET.append(waterRange-data[idx][1])
        
    WETDepth = []
    for idx, value in enumerate(WET):
        if (idx == 0):
            WETDepth.append(WET[idx]/2)
        else:
            WETDepth.append((WET[idx-1]+WET[idx])/2) 

    print("Water equivalent thicknesses")
    for idx, value in enumerate(WET):
        if idx > 0:
            print(f"CH: {idx+1}; WET: {WET[idx]-WET[idx-1]}")
        else:
            print(f"CH: {idx+1}; WET: {WET[idx]}")
    
    print(f"{WET}")
    print("Water equivalent Depth")
    print(WETDepth)


targetColorMap = ["#000000","#1f77b4", "#4e79a7", "#76b7b2", "#bab0ac", "#f28e2b", "#e15759", "#9c755f"]

path = "pbwo4" #h2o pbwo4 air DSB EJ256
name = "PbWO4" #H2O PbWO4 AIR DSB EJ-256
icru_el = 0 # 0=h20, 1=air 

colors = {
    "PbWO4":    (targetColorMap[0], targetColorMap[0]),
    "PTFE":  (targetColorMap[1], targetColorMap[1]),
    "H2O":   (targetColorMap[2], targetColorMap[2]),
    "Al": (targetColorMap[3], targetColorMap[3]),    # points, fit
}

def main():
    data = []
    dataWater = []

    layers = range(1,33)
    energies = [221, 225, 230, 235, 240, 245]

    energies = [180, 185, 190, 195, 200, 245]



    for index, value in enumerate(layers):
        if index == 18:
            data.append(get_range_energy(value, "pbwo4target", enable_output=True, enable_plot=True))
        else:
            data.append(get_range_energy(value, "pbwo4target", enable_output=True, enable_plot=False))
            
    for index, value in enumerate(energies):        
        dataWater.append(get_range_energy(value, "h2otarget", enable_output=True, enable_plot=False))
        
    # GetWET(data, dataWater)
    GetWETtarget(data, dataWater)

if __name__ == "__main__":
    main()