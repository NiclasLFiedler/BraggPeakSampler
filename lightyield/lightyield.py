import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Qt5Agg')

class BPSData:
    def __init__(self):
        self.channel = []
        self.count = []

    def add_data(self, channel, count):
        self.channel.append(channel)
        self.count.append(count)
    
    def clear(self):
        self.channel = []
        self.count = []

def process_files():
    bps_list = []

    # Loop over 15 files named BPS_CH$(number).txt
    path = "lateral"
    for i in range(0, 38):
        filename = f'data/2504/{path}/bps_pwo_lateral_{i}_na22.his.txt'
        if(i == 36):
            filename = f'data/2504/{path}/bps_pwo_lateral_window_0_na22.his.txt'
        elif(i == 37):
            filename = f'data/2504/{path}/bps_pwo_lateral_window_21_na22.his.txt'
        bps_data = BPSData()  # Create a new BPSData object for each file
        if os.path.exists(filename):
            with open(filename, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if line.strip():  # Check if line is not empty
                        channel, count = map(int, line.split())
                        bps_data.add_data(channel, count)
        else:
            print(f"File {filename} does not exist.")
        
        bps_list.append(bps_data)
        
    return bps_list
    
def gaussian(x, a, mu, sigma):
    return a * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Function to fit second peak
def fit_second_peak(BPSData, plot=True):
    # Find all peaks
    channel = np.array(BPSData.channel)
    count = np.array(BPSData.count)
    peaks, _ = find_peaks(count, width=50)

    sorted_peaks = sorted(peaks)
    second_peak_idx = sorted_peaks[0]
    peak_channel = channel[second_peak_idx]

    window = 75
    window_right = 100
    mask = (channel > peak_channel - window) & (channel < peak_channel + window+window_right)
    x_fit = channel[mask]
    y_fit = count[mask]

    p0 = [np.max(y_fit), peak_channel, 5]

    try:
        popt, _ = curve_fit(gaussian, x_fit, y_fit, p0=p0)
    except RuntimeError:
        raise RuntimeError("Gaussian fit did not converge.")

    if plot:
        plt.figure(figsize=(8, 5))
        plt.step(channel, count, label="Data")
        plt.plot(x_fit, gaussian(x_fit, *popt), 'r--', label="Gaussian Fit")
        plt.axvline(popt[1], color='g', linestyle=':', label=f"Mean = {popt[1]:.2f}")
        plt.legend()
        plt.xlabel("Channel")
        plt.ylabel("Counts")
        plt.title("Second Peak Gaussian Fit")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return popt  # [amplitude, mean, stddev]

def plot_bps_data(bps_data, title):
    plt.figure(figsize=(10, 5))
    plt.plot(bps_data.channel, bps_data.count, linestyle='-')
    plt.title(title)
    plt.xlabel('Channel')
    plt.ylabel('Count')
    plt.grid(True)
    plt.show()

def print_fitted_parameters(popt):
    amplitude, center, width = popt
    print(f"Fitted Gaussian Parameters:")
    print(f"  Amplitude: {amplitude:.2f}")
    print(f"  Center:    {center:.2f}")
    print(f"  Width:     {width:.2f}")

def main():
    bps_list = process_files()
    LY = []
    LY_err = []
    for i, bps_data in enumerate(bps_list):
        popt = fit_second_peak(bps_data, False)
        print(f"File {i}:")
        print(f"{i}, {popt[0]}, {popt[1]}, {popt[2]}")
        LY.append(popt[1])
        LY_err.append(popt[2])

    print([float(x) for x in LY])
    print()
    print([float(x) for x in LY_err])

if __name__ == "__main__":
    main()
