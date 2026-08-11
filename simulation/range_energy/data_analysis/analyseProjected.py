import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uproot
from pathlib import Path


def gaussian(x, amplitude, mean, sigma):
    """Gaussian function."""
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def load_and_extract_ranges(filename):
    """
    Load ROOT file and extract projected ranges (last z-position of trackid 0).
    
    Returns array of final positions for each event.
    """
    with uproot.open(filename) as file:
        tree = file["braggsampler"]
        
        # Load data
        pos = tree["pos"].array(library="np")
        trackid = tree["trackid"].array(library="np")
        event = tree["event"].array(library="np")
        
        # Extract last position for trackid 0 in each event
        ranges = []
        current_event = -1
        last_pos = None
        
        for i in range(len(trackid)):
            if trackid[i] == 1:
                if event[i] != current_event:
                    # New event, save previous range
                    if last_pos is not None:
                        ranges.append(last_pos)
                    current_event = event[i]
                
                last_pos = pos[i]
        
        # Don't forget the last one
        if last_pos is not None:
            ranges.append(last_pos)
    
    return np.array(ranges)


def fit_and_plot(ranges, energy, output_dir="./results"):
    """
    Fit Gaussian to ranges and create plot.
    
    Returns dict with range, error, and sigma.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Histogram
    counts, bin_edges = np.histogram(ranges, bins=500)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Fit
    peak_idx = np.argmax(counts)
    popt, pcov = curve_fit(
        gaussian,
        bin_centers,
        counts,
        p0=[counts[peak_idx], bin_centers[peak_idx], np.std(ranges) / 2],
        maxfev=5000
    )
    
    amplitude, mean, sigma = popt
    mean_err, sigma_err = np.sqrt(np.diag(pcov))[1:3]
    
    # Plot
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.bar(bin_centers, counts, width=bin_centers[1]-bin_centers[0], 
           alpha=0.7, edgecolor='black', label='Data')
    
    x_fit = np.linspace(bin_centers.min(), bin_centers.max(), 500)
    ax.plot(x_fit, gaussian(x_fit, *popt), 'r-', linewidth=2, label='Gaussian Fit')
    ax.axvline(mean, color='green', linestyle='--', linewidth=2, 
               label=f'Range = {mean:.2f} ± {mean_err:.3f} mm')
    
    ax.set_xlabel('Position (mm)', fontsize=12)
    ax.set_ylabel('Counts', fontsize=12)
    ax.set_title(f'Projected Range - {energy} MeV', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.savefig(output_dir / f"range_{energy}MeV.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"{energy} MeV: Range = {mean:.3f} ± {mean_err:.3f} mm, σ = {sigma:.3f} ± {sigma_err:.3f} mm")
    
    return {
        'energy': energy,
        'range': mean,
        'range_err': mean_err,
        'sigma': sigma,
        'sigma_err': sigma_err
    }


def main():
    # ===== EDIT THIS LIST WITH YOUR FILES =====
    files_and_energies = [
        ("h2oproj/raw_data_50.root", 50),
        ("h2oproj/raw_data_150.root", 150),
        ("h2oproj/raw_data_200.root", 200),
        ("h2oproj/raw_data_220.root", 220),
    ]
    # ==========================================
    
    output_dir = "./h2oproj"
    results = []
    
    for filename, energy in files_and_energies:
        if not Path(filename).exists():
            print(f"WARNING: {filename} not found, skipping")
            continue
        
        print(f"Processing {filename}...")
        ranges = load_and_extract_ranges(filename)
        print(f"  Found {len(ranges)} events")
        
        result = fit_and_plot(ranges, energy, output_dir)
        results.append(result)
    
    # Save results
    energies = np.array([r['energy'] for r in results])
    ranges = np.array([r['range'] for r in results])
    range_errors = np.array([r['range_err'] for r in results])
    sigmas = np.array([r['sigma'] for r in results])
    sigma_errors = np.array([r['sigma_err'] for r in results])
    
    output_file = Path(output_dir) / "results.npz"
    np.savez(output_file, 
             energies=energies, 
             ranges=ranges, 
             range_errors=range_errors,
             sigmas=sigmas,
             sigma_errors=sigma_errors)
    
    print(f"\nSaved to {output_file}")
    
    # Summary plot
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.errorbar(energies, ranges, yerr=range_errors, fmt='o-', markersize=8,
               capsize=5, linewidth=2)
    ax.set_xlabel('Beam Energy (MeV)', fontsize=12)
    ax.set_ylabel('Projected Range (mm)', fontsize=12)
    ax.set_title('Proton Projected Range vs Energy', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.savefig(Path(output_dir) / "range_vs_energy.pdf", bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main()