import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy import stats
import sys

def gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

def range_energy_relationship(energy, alpha, p):
    return alpha * energy**p


def analyze_geant4(df, layer_thickness=LAYER_THICKNESS):
    """
    Analyze GEANT4 data by binning hits by depth
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame from ROOT file
    layer_thickness : float
        Depth bin width in cm
    """
    
    print("\n=== GEANT4 Data Analysis ===")
    print(f"Total hits: {len(df)}")
    
    # Identify column names (handle various naming conventions)
    columns_lower = {col.lower(): col for col in df.columns}
    
    # Find depth column
    depth_col = None
    for candidate in ['z', 'pos.z', 'position.z', 'depth', 'posz']:
        if candidate in columns_lower:
            depth_col = columns_lower[candidate]
            break
    
    if not depth_col:
        print(f"Could not find depth column. Available columns:")
        for col in df.columns:
            print(f"  - {col}")
        return None
    
    # Find angle column
    angle_col = None
    for candidate in ['theta', 'angle', 'scattering_angle']:
        if candidate in columns_lower:
            angle_col = columns_lower[candidate]
            break
    
    if not angle_col:
        print(f"Could not find angle column. Available columns: {list(df.columns)}")
        return None
    
    # Find energy column
    energy_col = None
    for candidate in ['ekin', 'energy', 'kinetic_energy', 'e_kin']:
        if candidate in columns_lower:
            energy_col = columns_lower[candidate]
            break
    
    print(f"\nIdentified columns:")
    print(f"  Depth: {depth_col}")
    print(f"  Angle: {angle_col}")
    print(f"  Energy: {energy_col}")
    
    # Create depth bins
    depth_min = df[depth_col].min()
    depth_max = df[depth_col].max()
    depth_bins = np.arange(depth_min, depth_max + layer_thickness, layer_thickness)
    
    print(f"\nDepth range: {depth_min:.4f} - {depth_max:.4f} cm")
    print(f"Number of bins: {len(depth_bins)-1}")
    
    # Bin data by depth
    df['depth_bin'] = pd.cut(df[depth_col], bins=depth_bins)
    
    # Calculate statistics per depth bin
    geant4_results = []
    
    for bin_label, group in df.groupby('depth_bin', observed=True):
        if len(group) == 0:
            continue
        
        depth_center = bin_label.mid
        
        # Get angles for this bin
        angles = group[angle_col].values
        
        # Cumulative scattering angle (quadrature sum)
        cumulative_angle = np.sqrt(np.sum(angles**2))
        
        # Statistics
        mean_angle = np.mean(angles)
        std_angle = np.std(angles)
        n_hits = len(group)
        mean_energy = np.mean(group[energy_col].values) if energy_col else 0
        
        geant4_results.append({
            'depth_cm': depth_center,
            'cumulative_angle_mrad': cumulative_angle,
            'individual_angle_mrad': mean_angle,
            'std_angle_mrad': std_angle,
            'n_hits': n_hits,
            'mean_energy_MeV': mean_energy
        })
    
    geant4_df = pd.DataFrame(geant4_results)
    geant4_df = geant4_df.sort_values('depth_cm').reset_index(drop=True)
    
    print(f"Resulting bins with data: {len(geant4_df)}")
    
    return geant4_df


def fit_range_peak(ranges, bins=200):

    hist, edges = np.histogram(ranges, bins=bins)

    centers = 0.5*(edges[:-1]+edges[1:])
    peak_idx = np.argmax(hist)

    peak_pos = centers[peak_idx]


    # estimate sigma from distribution width
    sigma_guess = np.std(ranges)

    # fit window around peak
    fit_range = 3*sigma_guess

    mask = (centers > peak_pos-fit_range) & (centers < peak_pos+fit_range)


    p0 = [hist[peak_idx], peak_pos, sigma_guess]


    popt, pcov = curve_fit(gauss, centers[mask], hist[mask], p0=p0)
    perr = np.sqrt(np.diag(pcov))

    A, mu, sigma = popt

    return mu, sigma, A, centers, hist, popt, perr

materials = ["pbwo4proj", "h2oproj"]
# materials = ["pbwo4", "h2o"]

energies = [20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200]
energies = [220]
useMaterial = 1
material = materials[useMaterial]

files = [f"{material}/raw_data_{energy}.root" for energy in energies]

meanRanges = []

for file in files:
    tree = uproot.open(file)["braggsampler"]

    df = tree.arrays(
        ["event","pos","eDep", "theta", "trackid","eKin","eTot"],
        library="pd"
    )

    df = df[df.trackid == 1]
    
    depth_bin_width = 0.1  # mm

    max_depth = df.pos.max()

    depth_bins = np.arange(0, max_depth + depth_bin_width, depth_bin_width)

    theta_hist = []
    depth_centers = []
    
    theta_per_depth = [[] for _ in range(len(depth_bins)-1)]

    for event, g in df.groupby("event"):
        
        g = g.sort_values("pos")

        bins = np.digitize(g.pos.values, depth_bins)-1

        # keep only the last theta per depth bin
        unique_bins = np.unique(bins)

        for b in unique_bins:

            if 0 <= b < len(theta_per_depth):
                theta_per_depth[b].append(
                    g.theta.values[bins == b][-1]
                )


    theta_hist = []
    depth_bin_width = 0.1  # mm

    depth_bins = np.arange(
        0,
        df.pos.max()+depth_bin_width,
        depth_bin_width
    )

    # assign each step to a depth bin
    depth_index = np.digitize(df.pos.values, depth_bins)-1

    # collect RMS theta in each bin
    theta_hist = []
    depth_centers = []

    for i in range(len(depth_bins)-1):

        theta_values = df.theta.values[depth_index == i]

        if len(theta_values) > 0:
            theta_hist.append(
                np.sqrt(np.mean(theta_values**2))
            )
        else:
            theta_hist.append(np.nan)

        depth_centers.append(
            (depth_bins[i]+depth_bins[i+1])/2
        )

plt.figure(figsize=(8,5))

plt.plot(
    depth_centers,
    np.array(theta_hist)*1000
)

plt.xlabel("Depth [mm]")
plt.ylabel(r"$\theta_x$ RMS [mrad]")
plt.grid()

plt.show()