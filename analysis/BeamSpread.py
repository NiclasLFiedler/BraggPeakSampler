import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define Gaussian
def gaus(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))



def main():
    # Open ROOT file
    file = uproot.open("histograms.root")

    # Select histograms that start with "h_enty_"
    entry_hists = {key: file[key] for key in file.keys() if key.startswith("h_enty_")}

    # Start a single figure
    plt.figure(figsize=(10, 6))

    for name, hist in entry_hists.items():
        # Extract histogram data
        counts, bin_edges = hist.to_numpy()
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # Fit to Gaussian
        try:
            popt, _ = curve_fit(gaus, bin_centers, counts, p0=[np.max(counts), bin_centers[np.argmax(counts)], 1.0])
            A, mu, sigma = popt
            print(f"{name}: A={A:.2f}, μ={mu:.2f}, σ={sigma:.2f}")
        except RuntimeError:
            print(f"❌ Fit failed for {name}")
            continue

        # Plot histogram
        plt.step(bin_edges[:-1], counts, where="post", label=f"{name} hist")

        # Plot Gaussian fit
        x_fit = np.linspace(bin_centers[0], bin_centers[-1], 1000)
        plt.plot(x_fit, gaus(x_fit, *popt), '--', label=f"{name} fit")

    plt.xlabel("EntryX / cm")
    plt.ylabel("Counts")
    plt.title("All Histograms with Gaussian Fits")
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
else:
    print("This script is intended to be run as a standalone program.")