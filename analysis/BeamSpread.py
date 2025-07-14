import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator, FuncFormatter
# Define Gaussian
def gaus(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def format_pi(x, pos):
    frac = x / np.pi
    if np.isclose(frac, 0):
        return "0"
    elif np.isclose(frac, 1):
        return r"$\pi$"
    elif np.isclose(frac, -1):
        return r"$-\pi$"
    elif np.isclose(frac, 2):
        return r"$2\pi$"
    elif np.isclose(frac, -2):
        return r"$-2\pi$"
    else:
        # Reduce fraction to nice format
        from fractions import Fraction
        frac = Fraction(frac).limit_denominator(8)  # π/8 resolution
        num, denom = frac.numerator, frac.denominator
        if denom == 1:
            return rf"${num}\pi$"
        else:
            return rf"$\frac{{{num}\pi}}{{{denom}}}$"

def radians_to_degrees(x):
    return np.degrees(x)

def degrees_to_radians(x):
    return np.radians(x)

def main():
    # Open ROOT file
    file = uproot.open("histograms.root")

    # Select histograms that start with "h_enty_"
    entry_hists = {key: file[key] for key in file.keys() if key.startswith("h_entry_")}
    print(f"Found {len(entry_hists)} histograms starting with 'h_entry_'")

    # Start a single figure
    plt.figure(figsize=(10, 6))

    positions = [6.77195, 20.8666, 35.2395, 49.7721, 64.4184, 79.1533, 93.9607, 108.83, 123.753, 138.723, 153.736, 166.459, 176.875, 187.307, 197.755, 208.216, 218.691, 229.178, 239.678, 250.19, 260.713, 271.247, 281.791, 292.346, 302.91, 313.484, 324.067, 334.659, 345.259, 355.868, 366.485, 377.11]
    popts = []
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
            popts.append([0, 0, 0])
            continue
        popts.append(popt)
        # Plot histogram
        plt.step(bin_edges[:-1], counts, where="post", label=f"{name} hist")

        # Plot Gaussian fit
        x_fit = np.linspace(bin_centers[0], bin_centers[-1], 1000)
        plt.plot(x_fit, gaus(x_fit, *popt), '--', label=f"{name} fit")

    plt.xlabel("EntryX / cm")
    plt.ylabel("Counts")
    plt.title("All Histograms with Gaussian Fits")
    #plt.legend()
    
    plt.show()

    print("len popsts:", len(popts))
    plt.plot(positions, [popt[2] for popt in popts], 'o-')
    plt.xlabel("Position / cm")
    plt.ylabel("Std Dev of Gaussian Fit")
    plt.title("Std Dev Position of Gaussian Fits")
    plt.grid()
    plt.show()

    # Select histograms that start with "h_enty_"
    entry_hists = {key: file[key] for key in file.keys() if key.startswith("h_angle_")}
    fig, ax = plt.subplots()

    ax.xaxis.set_major_locator(MultipleLocator(base=np.pi/8))
    ax.xaxis.set_major_formatter(FuncFormatter(format_pi))
    
    positions = [6.87413, 20.8958, 35.2571, 49.7848, 64.4284, 79.1614, 93.9677, 108.836, 123.758, 138.728, 153.74, 166.46, 176.877, 187.309, 197.756, 208.217, 218.692, 229.18, 239.679, 250.191, 260.714, 271.248, 281.793, 292.347, 302.911, 313.485, 324.068, 334.66, 345.26, 355.869, 366.486, 377.111]    

    popts = []


    i = 0
    for name, hist in entry_hists.items():
        i += 1
        #if i > 2 :
        #     break
        # Extract histogram data
        counts, bin_edges = hist.to_numpy()
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        ax.step(bin_edges[:-1], counts, where="post", label=f"{name} hist")

    secax = ax.secondary_xaxis('top', functions=(radians_to_degrees, degrees_to_radians))
    secax.set_xlabel('Angle (degrees)')
    ax.set_xlabel('Angle (radians)')
    ax.set_ylabel("Counts")
    ax.grid(True)
    ax.set_xlim(0, np.pi) 
    ax.set_yscale('log')
    plt.show()

if __name__ == "__main__":
    main()
else:
    print("This script is intended to be run as a standalone program.")