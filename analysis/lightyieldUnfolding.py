# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
# ---

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import curve_fit

# %%
from scipy.stats import poisson
from scipy.signal import convolve
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import ROOT
import uproot
import matplotlib
matplotlib.use('TkAgg')  # or 'Agg' for non-interactive
import os
ROOT.gSystem.Load("libRooUnfold")
from concurrent.futures import ProcessPoolExecutor

def numpy_to_th1(name, title, values, edges):
    values = np.ascontiguousarray(values, dtype=np.float64)
    edges = np.ascontiguousarray(edges, dtype=np.float64)

    n_bins = len(values)

    h = ROOT.TH1D(name, title, n_bins, edges)
    for i, v in enumerate(values):
        h.SetBinContent(i+1, v)
        h.SetBinError(i+1, np.sqrt(v))  # Poisson errors
    return h

def numpy_to_th2(name, title, R, xedges, yedges):
    
    R = np.ascontiguousarray(R, dtype=np.float64)
    xedges = np.ascontiguousarray(xedges, dtype=np.float64)
    yedges = np.ascontiguousarray(yedges, dtype=np.float64)    
    
    h = ROOT.TH2D(
        name, title,
        len(xedges)-1, xedges,
        len(yedges)-1, yedges
    )
    for ix in range(R.shape[0]):      # measured
        for iy in range(R.shape[1]):  # true
            h.SetBinContent(ix+1, iy+1, R[ix, iy])
    return h

def calcBinEdges(centers):
    centers = np.array(centers)
    edges = np.zeros(len(centers) + 1)
    edges[1:-1] = (centers[1:] + centers[:-1]) / 2
    edges[0] = centers[0] - (centers[1] - centers[0]) / 2
    edges[-1] = centers[-1] + (centers[-1] - centers[-2]) / 2
    return edges

def th1_to_numpy(h):
    """Return bin centers and bin contents"""
    bins = np.array([h.GetBinCenter(i+1) for i in range(h.GetNbinsX())])
    contents = np.array([h.GetBinContent(i+1) for i in range(h.GetNbinsX())])
    return bins, contents

# %%
class BayesianUnfoldingWithDRF:
    """
    Bayesian Unfolding with a functional Detector Response Function (DRF)
    """
    
    def __init__(self, n_true_bins=4001, n_measured_bins=4001, n_iterations=10):
        self.n_true_bins = n_true_bins
        self.n_measured_bins = n_measured_bins
        self.n_iterations = n_iterations

    def drf_function(self, measured_energy, true_energy, resolution=0.15, tail_strength=0.2):
        """
        EXAMPLE DRF function - REPLACE THIS WITH YOUR ACTUAL DRF!
        
        This is a realistic DRF for a scintillator-SiPM system:
        - Gaussian main peak with energy-dependent resolution
        - Low-energy tail from incomplete light collection
        
        Parameters:
        measured_energy: the observed energy (MeV)
        true_energy: the true deposited energy (MeV)
        resolution: energy resolution (σ/E) at 1 MeV
        tail_strength: relative strength of low-energy tail
        
        Returns:
        probability density for measuring 'measured_energy' given 'true_energy'
        """
        # Energy resolution typically scales as 1/sqrt(E)
        sigma = resolution * np.sqrt(true_energy) * true_energy
        
        # Main Gaussian peak
        gaussian = np.exp(-0.5 * ((measured_energy - true_energy) / sigma)**2) / (sigma * np.sqrt(2 * np.pi))
        
        # Low-energy tail (common in scintillators due to light collection effects)
        if measured_energy < true_energy:
            tail = tail_strength * np.exp(-(true_energy - measured_energy) / (0.1 * true_energy)) / true_energy
        else:
            tail = 0
            
        return gaussian + tail
    
    def build_response_matrix(self, true_energy_range, measured_energy_range, custom_drf=None):
        """
        Build response matrix from DRF function
        
        Parameters:
        true_energy_range: tuple (min, max) for true energies
        measured_energy_range: tuple (min, max) for measured energies  
        custom_drf: your custom DRF function. If None, uses the example above
        
        Returns:
        response_matrix: R[i,j] = P(measured_bin i | true_bin j)
        true_energies: center of true energy bins
        measured_energies: center of measured energy bins
        """
        # Use custom DRF if provided, otherwise use example
        drf = custom_drf if custom_drf is not None else self.drf_function
        
        # Create energy bins
        self.true_edges = np.linspace(true_energy_range[0], true_energy_range[1], self.n_true_bins + 1)
        self.measured_edges = np.linspace(measured_energy_range[0], measured_energy_range[1], self.n_measured_bins + 1)
        
        self.true_centers = (self.true_edges[1:] + self.true_edges[:-1]) / 2
        self.measured_centers = (self.measured_edges[1:] + self.measured_edges[:-1]) / 2
        bin_widths = self.measured_edges[1:] - self.measured_edges[:-1]
        
        response_matrix = np.zeros((self.n_measured_bins, self.n_true_bins))
        
        print("Building response matrix...")
        for j, true_energy in enumerate(self.true_centers):
            if j % 100 == 0:
                print(f"  Processing true energy bin {j+1}/{self.n_true_bins}")
            
            # For each true energy, calculate probability in each measured bin
            # for i in range(self.n_measured_bins):
            #     # Integrate DRF over the measured energy bin
            #     low = self.measured_edges[i]
            #     high = self.measured_edges[i + 1]
                
            #     # Numerical integration of DRF over the bin
            #     integral, error = quad(drf, low, high, args=(true_energy))
            #     response_matrix[i, j] = integral
            
            energyResponse = drf(true_energy, self.measured_centers)
            
            response_matrix[:, j] = energyResponse * bin_widths
        
        # Normalize each column (true energy) to sum to 1
        response_matrix = response_matrix / np.sum(response_matrix, axis=0, keepdims=True)

        self.R = response_matrix
        return response_matrix, self.true_centers, self.measured_centers
    
    def load_response_matrix(self, true_energy_range, measured_energy_range, filename):
        self.true_edges = np.linspace(true_energy_range[0], true_energy_range[1], self.n_true_bins + 1)
        self.measured_edges = np.linspace(measured_energy_range[0], measured_energy_range[1], self.n_measured_bins + 1)
        
        self.true_centers = (self.true_edges[1:] + self.true_edges[:-1]) / 2
        self.measured_centers = (self.measured_edges[1:] + self.measured_edges[:-1]) / 2
        file = uproot.open(filename)
        hist = file["h_responseMatrix"]
        response_matrix, _, _ = hist.to_numpy()
        
        self.R = response_matrix
        return response_matrix, self.true_centers, self.measured_centers
        
    def unfold(self, measured_spectrum, prior=None, efficiency=None):
        """
        Unfold measured spectrum using Bayesian method
        
        Parameters:
        measured_spectrum: 1D array of measured counts (binned according to measured_edges)
        prior: initial guess for true spectrum
        efficiency: detection efficiency for each true bin
        
        Returns:
        unfolded_spectrum: estimated true spectrum
        """
        if not hasattr(self, 'R'):
            raise ValueError("Must build response matrix first using build_response_matrix()")
        
        g = measured_spectrum.astype(float)
        f = prior.copy() if prior is not None else np.ones(self.n_true_bins)
        
        # Normalize prior
        f = f / np.sum(f) * np.sum(g)
        
        if efficiency is None:
            efficiency = np.sum(self.R, axis=0)
        
        max_index = np.argmax(measured_spectrum)

        self.history = [f.copy()]
        self.chi2 = [0]
        self.closure_errors = [0]
        print("\nStarting Bayesian unfolding...")
        for iteration in range(self.n_iterations):
            g_expected = np.dot(self.R, f)
            
            g_expected[g_expected == 0] = 1e-10
            
            correction_factor = np.dot(self.R.T, g / g_expected)
            f_new = f * correction_factor / efficiency
            f_unsmoothed = f_new.copy()
            f_new = self._smooth_spectrum(f_new, strength=0.1)
            # f_new = self._tikhonov_regularize(f_new, f, 0.01)
            f_new = np.maximum(f_new, 0)
            
            # mask = (self.history[0] > 0) & (f_new > 0)
            # chi2 = np.sum((f_new[mask] - self.history[0][mask])**2 / self.history[0][mask])

            predicted_measured = np.dot(self.R, f_new)
            mask = (measured_spectrum > 0) & (predicted_measured > 0)
            mask[:max_index+130] = False
            n_dof = np.sum(mask) - 1
            chi2 = np.sum((predicted_measured[mask] - measured_spectrum[mask])**2 / measured_spectrum[mask])/n_dof
            closure_error = np.sqrt(np.mean(((measured_spectrum[mask] - predicted_measured[mask]) / measured_spectrum[mask])**2))
            self.chi2.append(chi2)
            self.closure_errors.append(closure_error)
            
            if(self.chi2[-2]-self.chi2[-1]<0.03 and iteration>3):
                print(f"Converged after {iteration+1} iterations")
                f = f_unsmoothed
                self.history.append(f.copy())
                
                break

            f = f_new
            self.history.append(f.copy())
                                
        self.unfolded = f
        return f

    def unfold_regularized(self, measured_spectrum, prior=None, regularization_strength=0.01, early_stopping=True):
        """
        Unfold with Tikhonov regularization for better convergence
        """
        if not hasattr(self, 'R'):
            raise ValueError("Must build response matrix first")

        g = measured_spectrum.astype(float)
        f = prior.copy() if prior is not None else np.ones(self.n_true_bins)

        # Normalize prior
        f = f / np.sum(f) * np.sum(g)

        efficiency = np.sum(self.R, axis=0)

        self.history = [f.copy()]
        self.closure_errors = []

        print("\nStarting regularized Bayesian unfolding...")
        for iteration in range(self.n_iterations):
            g_expected = np.dot(self.R, f)
            g_expected[g_expected == 0] = 1e-10

            correction_factor = np.dot(self.R.T, g / g_expected)
            f_new = f * correction_factor / efficiency

            f_new = self._tikhonov_regularize(f_new, f, regularization_strength)
            f_new = np.maximum(f_new, 0)

            relative_change = np.sum(np.abs(f_new - f)) / np.sum(f)

            f = f_new
            self.history.append(f.copy())

            pred_meas = np.dot(self.R, f)
            closure_err = np.sqrt(np.mean(((g - pred_meas) / (g + 1e-10))**2))
            self.closure_errors.append(closure_err)

            print(f"Iteration {iteration+1}: Total counts = {np.sum(f):.1f}, "
                  f"Change = {relative_change:.2e}, Closure = {closure_err:.4f}")

            # Early stopping
            if early_stopping and relative_change < 1e-4 and iteration > 3:
                print(f"Converged after {iteration+1} iterations")
                break
            
        self.unfolded = f
        return f

    def _tikhonov_regularize(self, f_new, f_old, strength):
        """Apply Tikhonov regularization to reduce oscillations"""
        # Simple gradient penalty: minimize curvature
        from scipy import sparse

        # Create second derivative matrix
        n = len(f_new)
        D2 = sparse.diags([1, -2, 1], [0, 1, 2], shape=(n-2, n)).toarray()

        # Calculate regularization term
        reg_term = strength * np.dot(D2.T, np.dot(D2, f_new))

        # Apply regularization (subtract gradient)
        f_regularized = f_new - reg_term

        return np.maximum(f_regularized, 0)

    def _smooth_spectrum(self, spectrum, strength=0.1):
        """Apply mild smoothing to reduce oscillations"""
        from scipy.ndimage import gaussian_filter1d
        return (1 - strength) * spectrum + strength * gaussian_filter1d(spectrum, sigma=10.0)
    
    def calculate_closure(self, measured_spectrum):
        """Calculate how well the unfolded spectrum reproduces the measurement"""
        if hasattr(self, 'unfolded'):
            predicted_measured = np.dot(self.R, self.unfolded)
            
            # Avoid division by zero
            mask = measured_spectrum > 0
            
            # Calculate normalized chi-squared
            chi2 = np.sum((measured_spectrum[mask] - predicted_measured[mask])**2 / measured_spectrum[mask])
            n_dof = np.sum(mask) - 1  # degrees of freedom
            
            # Calculate relative closure error
            closure_error = np.sqrt(np.mean(((measured_spectrum[mask] - predicted_measured[mask]) / measured_spectrum[mask])**2))
            
            return closure_error, chi2/n_dof, predicted_measured
        else:
            raise ValueError("Must run unfold() first")
    
    def calculate_closure_fixed(self, measured_spectrum):
        """Proper closure calculation that handles zero bins correctly"""
        if hasattr(self, 'unfolded'):
            predicted_measured = np.dot(self.R, self.unfolded)

            # Method 1: Use only bins with sufficient counts
            threshold = np.max(measured_spectrum) * 0.01  # 1% of peak
            mask = measured_spectrum > threshold

            if np.sum(mask) < 5:  # If too few bins, use non-zero bins
                mask = measured_spectrum > 0

            # Calculate metrics only on meaningful bins
            if np.sum(mask) > 0:
                # Relative RMS error on meaningful bins
                rel_errors = np.abs(measured_spectrum[mask] - predicted_measured[mask]) / measured_spectrum[mask]
                closure_rms = np.sqrt(np.mean(rel_errors**2))

                # Chi-squared per DOF
                chi2 = np.sum((measured_spectrum[mask] - predicted_measured[mask])**2 / measured_spectrum[mask])
                chi2_per_dof = chi2 / len(mask)

                # Correlation
                correlation = np.corrcoef(measured_spectrum[mask], predicted_measured[mask])[0, 1]
            else:
                closure_rms = float('inf')
                chi2_per_dof = float('inf')
                correlation = 0

            # Method 2: Overall shape metrics (less sensitive to zeros)
            total_measured = np.sum(measured_spectrum)
            total_predicted = np.sum(predicted_measured)
            total_error = np.abs(total_measured - total_predicted) / total_measured

            # Method 3: Wasserstein distance (shape similarity)
            from scipy.stats import wasserstein_distance
            wasserstein = wasserstein_distance(measured_spectrum, predicted_measured)

            print(f"\n=== ROBUST CLOSURE METRICS ===")
            print(f"Bins used for metrics: {np.sum(mask)}/{len(measured_spectrum)}")
            print(f"Relative RMS Error: {closure_rms:.4f} ({closure_rms*100:.1f}%)")
            print(f"Chi-squared per DOF: {chi2_per_dof:.2f}")
            print(f"Correlation: {correlation:.4f}")
            print(f"Total Counts Error: {total_error:.4f} ({total_error*100:.1f}%)")
            print(f"Wasserstein Distance: {wasserstein:.4f}")

            return {
                'rms_error': closure_rms,
                'chi2_per_dof': chi2_per_dof,
                'correlation': correlation,
                'total_error': total_error,
                'wasserstein': wasserstein,
                'predicted': predicted_measured,
                'mask': mask
            }

    def check_response_matrix(self):
        """Diagnose potential issues with the response matrix"""
        print("\nResponse Matrix Diagnostics:")
        print(f"Shape: {self.R.shape}")
        print(f"Min value: {np.min(self.R):.2e}")
        print(f"Max value: {np.max(self.R):.2e}")
        print(f"Matrix condition number: {np.linalg.cond(self.R):.2e}")

        # Check for columns with very small sums (ill-conditioned)
        column_sums = np.sum(self.R, axis=0)
        zero_columns = np.sum(column_sums < 1e-10)
        print(f"Columns with near-zero sum: {zero_columns}")

        if zero_columns > 0:
            print("WARNING: Response matrix has columns with near-zero sum!")
            print("This can cause convergence issues.")

    def calculate_closure_robust(self, measured_spectrum):
        """Robust closure calculation that handles edge cases"""
        if hasattr(self, 'unfolded'):
            predicted_measured = np.dot(self.R, self.unfolded)

            # Method 1: Normalized RMS (handles zeros)
            diff = measured_spectrum - predicted_measured
            rms_error = np.sqrt(np.mean(diff**2))
            normalization = np.sqrt(np.mean(measured_spectrum**2))
            closure_rms = rms_error / normalization if normalization > 0 else float('inf')

            # Method 2: Relative error in total counts
            total_measured = np.sum(measured_spectrum)
            total_predicted = np.sum(predicted_measured)
            closure_total = np.abs(total_measured - total_predicted) / total_measured

            # Method 3: Correlation coefficient
            correlation = np.corrcoef(measured_spectrum, predicted_measured)[0, 1]

            # Method 4: Exclude very low count bins
            mask = measured_spectrum > np.max(measured_spectrum) * 0.01  # Exclude lowest 1%
            if np.sum(mask) > 10:  # Only use if we have enough bins
                rel_error_masked = np.sqrt(np.mean(((measured_spectrum[mask] - predicted_measured[mask]) / measured_spectrum[mask])**2))
            else:
                rel_error_masked = float('inf')

            print(f"\nClosure Diagnostics:")
            print(f"RMS Error: {closure_rms:.4f} ({closure_rms*100:.1f}%)")
            print(f"Total Counts Error: {closure_total:.4f} ({closure_total*100:.1f}%)")
            print(f"Correlation: {correlation:.4f}")
            print(f"Masked Relative Error: {rel_error_masked:.4f} ({rel_error_masked*100:.1f}%)")

            return {
                'rms_error': closure_rms,
                'total_error': closure_total, 
                'correlation': correlation,
                'masked_error': rel_error_masked,
                'predicted': predicted_measured
            }
        else:
            raise ValueError("Must run unfold() first")

    def debug_closure_issue(self, measured_spectrum, predicted_measured):
        """Debug why closure error is huge but plots look good"""
        print("\n=== CLOSURE DEBUGGING ===")
        print(f"Measured spectrum sum: {np.sum(measured_spectrum):.1f}")
        print(f"Predicted spectrum sum: {np.sum(predicted_measured):.1f}")
        print(f"Ratio (predicted/measured): {np.sum(predicted_measured)/np.sum(measured_spectrum):.3f}")

        # Check individual bins
        diff = measured_spectrum - predicted_measured
        rel_error = np.abs(diff) / (measured_spectrum + 1e-10)  # Avoid division by zero

        print(f"\nTop 5 largest relative errors:")
        worst_indices = np.argsort(rel_error)[-5:][::-1]
        for idx in worst_indices:
            print(f"  Bin {idx}: Measured={measured_spectrum[idx]:.1f}, "
                  f"Predicted={predicted_measured[idx]:.1f}, "
                  f"Error={rel_error[idx]*100:.1f}%")

        # Check for zeros in measured spectrum
        zero_bins = np.sum(measured_spectrum == 0)
        print(f"\nBins with zero measured counts: {zero_bins}")

        if zero_bins > 0:
            print("WARNING: Zero bins in measured spectrum can inflate closure error!")


# %%
def save_response_matrix_th2d(response_matrix, true_edges, measured_edges, 
                             filename="response_matrix.root", histname="response_matrix"):
    """
    Save response matrix as TH2D in ROOT file
    
    Parameters:
    response_matrix: 2D numpy array (n_measured_bins, n_true_bins)
    true_edges: 1D array of true energy bin edges
    measured_edges: 1D array of measured energy bin edges  
    filename: output ROOT filename
    histname: name of the TH2D histogram
    """
    
    # Convert edges to ROOT arrays
    true_edges_root = np.array('d', true_edges)
    measured_edges_root = np.array('d', measured_edges)
    
    n_true_bins = len(true_edges) - 1
    n_measured_bins = len(measured_edges) - 1
    
    # Create TH2D
    h_response = ROOT.TH2D(histname, "Detector Response Matrix",
                          n_true_bins, true_edges_root,
                          n_measured_bins, measured_edges_root)
    
    # Set axis titles
    h_response.GetXaxis().SetTitle("True Energy")
    h_response.GetYaxis().SetTitle("Measured Energy")
    h_response.GetZaxis().SetTitle("Probability")
    
    # Fill the histogram
    for i in range(n_measured_bins):      # y-axis (measured)
        for j in range(n_true_bins):      # x-axis (true)
            # TH2D bin indexing: (x_bin, y_bin) where x=true, y=measured
            h_response.SetBinContent(j+1, i+1, response_matrix[i, j])
    
    # Save to file
    output_file = ROOT.TFile(filename, "RECREATE")
    h_response.Write()
    output_file.Close()
    
    print(f"Response matrix saved to {filename}")
    print(f"Histogram name: {histname}")
    print(f"Dimensions: {n_true_bins} true bins × {n_measured_bins} measured bins")
    
    return h_response

# %%
def save_comprehensive_response_data(response_matrix, true_edges, measured_edges, 
                                   true_centers, measured_centers, efficiency,
                                   filename="response_data.root"):
    """
    Save comprehensive response matrix data including 1D projections
    """
    # Convert edges to ROOT arrays
    true_edges_root = np.array('d', true_edges)
    measured_edges_root = np.array('d', measured_edges)
    true_centers_root = np.array('d', true_centers)
    measured_centers_root = np.array('d', measured_centers)
    
    n_true_bins = len(true_edges) - 1
    n_measured_bins = len(measured_edges) - 1
    
    # Create output file
    output_file = ROOT.TFile(filename, "RECREATE")
    
    # 1. Main response matrix (TH2D)
    h_response = ROOT.TH2D("response_matrix", "Detector Response Matrix R(i,j)=P(measured i|true j)",
                          n_true_bins, true_edges_root,
                          n_measured_bins, measured_edges_root)
    h_response.GetXaxis().SetTitle("True Energy [MeV]")
    h_response.GetYaxis().SetTitle("Measured Energy [MeV]")
    h_response.GetZaxis().SetTitle("Probability Density")
    
    # 2. Efficiency curve (TH1D)
    h_efficiency = ROOT.TH1D("efficiency", "Detection Efficiency",
                            n_true_bins, true_edges_root)
    h_efficiency.GetXaxis().SetTitle("True Energy [MeV]")
    h_efficiency.GetYaxis().SetTitle("Efficiency")
    
    # 3. Example response slices (TH1D for different true energies)
    h_response_slices = []
    
    # Fill histograms
    for i in range(n_measured_bins):
        for j in range(n_true_bins):
            h_response.SetBinContent(j+1, i+1, response_matrix[i, j])
    
    for j in range(n_true_bins):
        h_efficiency.SetBinContent(j+1, efficiency[j])
        h_efficiency.SetBinError(j+1, 0.01)  # Small error for plotting
    
    # Create example slices at different true energies
    slice_positions = [0.25, 0.5, 0.75]  # Fractions of energy range
    for frac in slice_positions:
        j = int(n_true_bins * frac)
        if j >= n_true_bins:
            j = n_true_bins - 1
        
        slice_name = f"response_slice_true_{true_centers[j]:.2f}MeV"
        h_slice = ROOT.TH1D(slice_name, f"Response at {true_centers[j]:.2f} MeV",
                           n_measured_bins, measured_edges_root)
        h_slice.GetXaxis().SetTitle("Measured Energy [MeV]")
        h_slice.GetYaxis().SetTitle("Probability Density")
        
        for i in range(n_measured_bins):
            h_slice.SetBinContent(i+1, response_matrix[i, j])
        
        h_slice.Write()
        h_response_slices.append(h_slice)
    
    # Write all histograms
    h_response.Write()
    h_efficiency.Write()
    
    # Create a directory for projections
    output_file.mkdir("projections")
    output_file.cd("projections")
    
    # Save projections
    h_projX = h_response.ProjectionX("projection_true")
    h_projY = h_response.ProjectionY("projection_measured")
    h_projX.Write()
    h_projY.Write()
    
    output_file.cd()  # Back to main directory
    
    # Save metadata as TTree
    tree = ROOT.TTree("metadata", "Response Matrix Metadata")
    
    # Create variables for the tree
    n_true = np.array([n_true_bins], dtype='i')
    n_measured = np.array([n_measured_bins], dtype='i')
    true_min = np.array([true_edges[0]], dtype='f')
    true_max = np.array([true_edges[-1]], dtype='f')
    measured_min = np.array([measured_edges[0]], dtype='f')
    measured_max = np.array([measured_edges[-1]], dtype='f')
    
    # Create branches
    tree.Branch("n_true_bins", n_true, "n_true_bins/I")
    tree.Branch("n_measured_bins", n_measured, "n_measured_bins/I")
    tree.Branch("true_energy_min", true_min, "true_energy_min/F")
    tree.Branch("true_energy_max", true_max, "true_energy_max/F")
    tree.Branch("measured_energy_min", measured_min, "measured_energy_min/F")
    tree.Branch("measured_energy_max", measured_max, "measured_energy_max/F")
    
    tree.Fill()
    tree.Write()
    
    output_file.Close()
    
    print(f"Comprehensive response data saved to {filename}")
    print("Contents:")
    print("  - response_matrix: TH2D main response matrix")
    print("  - efficiency: TH1D detection efficiency")
    print("  - response_slice_*: TH1D example response slices")
    print("  - projections/: directory with 1D projections")
    print("  - metadata: TTree with binning information")

# %%
def save_response_matrix(self, filename="response_matrix.root"):
    """Save the current response matrix to ROOT file"""
    if not hasattr(self, 'R'):
        print("Error: Build response matrix first using build_response_matrix()")
        return
    
    save_comprehensive_response_data(
        self.R, self.true_edges, self.measured_edges,
        self.true_centers, self.measured_centers,
        np.sum(self.R, axis=0),  # efficiency
        filename
    )

# %%
# Example usage with your data
if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(42)
    
    print("Bayesian Unfolding with Functional DRF")
    print("=" * 50)

    
    def gaussian(x, amplitude, mean, sigma):
        return amplitude * np.exp(-0.5 * ((x - mean) / sigma)**2)

    def pmt_response(true_pe, pe_axis, n_max=1000):
        """
        Calculate PMT response function for a given energy deposition
        Parameters:
        - true_adc: True adc value
        - mean_SPE: mean of single photoelectron peak
        - sigma_SPE: Standard deviation of single photoelectron peak
        - mean_ped: mean of pedastal
        - sigma_ped: Standard deviation of pedastal
        - pe_range: Tuple (min, max, n_points) for pe axis
        - n_max: Maximum number of photoelectrons to consider
        Returns:
        - pe_axis: Array of pe values
        - response: Probability density function
        - components: Individual photoelectron peaks for plotting
        """
        mean_SPE = 110.342
        sigma_SPE= 10.0004
        mean_ped= 85.5295
        sigma_ped= 0.48595
        PDE = 0.2316
        sigma_SPE = (sigma_SPE) / (mean_SPE-mean_ped)#/0.2316
        sigma_ped = (sigma_ped) / (mean_SPE-mean_ped)#/0.2316
        response = np.zeros_like(pe_axis)
        bin_width = unfolder.true_centers[1]-unfolder.true_centers[0]
        underflow_axis = np.arange(-(bin_width)*400+unfolder.true_centers[0], -unfolder.true_centers[0], unfolder.true_centers[1]-unfolder.true_centers[0])
        overflow_axis = np.arange(unfolder.true_centers[-1]+unfolder.true_centers[0], (bin_width)*400+unfolder.true_centers[-1], unfolder.true_centers[1]-unfolder.true_centers[0])

        underflow_response = np.zeros_like(underflow_axis)
        overflow_response = np.zeros_like(overflow_axis)
        components = []

        for n in range(0, n_max + 1):    
            poisson_prob = poisson.pmf(n, true_pe*PDE)
            if poisson_prob > 1e-10:  
                sigma_n = np.sqrt((np.sqrt(n) * sigma_SPE)**2+sigma_ped**2) if n > 0 else np.sqrt(sigma_SPE**2+sigma_ped**2)
                if sigma_n < 1e-10:
                    sigma_n = sigma_SPE
                gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                          np.exp(-0.5 * ((pe_axis*PDE - n) / sigma_n) ** 2)
                component = poisson_prob * gaussian

                if(component[0]>1e-5):
                    underflow_gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                              np.exp(-0.5 * ((underflow_axis*PDE - n) / sigma_n) ** 2)
                    underflow_component = poisson_prob * underflow_gaussian

                    underflow_response += underflow_component
                    component[0] += np.sum(underflow_component)
                
                if(component[-1]>1e-5):
                    overflow_gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                              np.exp(-0.5 * ((overflow_axis*PDE - n) / sigma_n) ** 2)
                    overflow_component = poisson_prob * overflow_gaussian

                    overflow_response += overflow_component
                    component[-1] += np.sum(overflow_component)
                components.append((n, n, sigma_n, component))
                response += component
                
        return response

    def findPeaks(unfolded_spectrum, measured_spectrum, true_energies, fheight=50, fdistance=30, fwidth=30):
        measuredPeaks, _ = find_peaks(measured_spectrum, height=1000, distance=30,width=3)
        print(measuredPeaks)
        print()
        peaks, _ = find_peaks(unfolded_spectrum, height=fheight, distance=fdistance, width=fwidth)
        print(peaks)
        if(len(measuredPeaks) > 0):
            peaks = peaks[peaks > (measuredPeaks[0]+100)]
        print(peaks)
        print()
        selected_peaks = peaks

        params = []
        x_fits = []
        x_fit = []
        popt = []
        for idx, p in enumerate(selected_peaks):
            if idx == 1:
                fwidth = 50
            mask = (true_energies > true_energies[p] - (fwidth-15)) & (true_energies < true_energies[p] + fwidth-15)
            x_fit = true_energies[mask]
            y_fit = unfolded_spectrum[mask]

            amp0 = unfolded_spectrum[p]
            mean0 = true_energies[p]
            sigma0 = 2.0

            popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[amp0, mean0, sigma0])
            perr = np.sqrt(np.diag(pcov))
            params.append([popt, perr])
            x_fits.append(x_fit)
        return selected_peaks, params, x_fits
    
    def saveGif():
        import imageio.v2 as imageio
        os.makedirs("frames", exist_ok=True)
        filenames = []
        for i, E in enumerate(true_energies):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(R[:,i], color='C0')
            ax.set_title(f"Response for Energy = {E:.1f} MeV")
            ax.set_xlabel("Detector Channel")
            ax.set_ylabel("Signal (a.u.)")
            ax.set_ylim(0, R[:,i].max()*1.1)

            fname = f"gif/frame_{i:03d}.png"
            plt.savefig(fname)
            plt.close(fig)
            filenames.append(fname)

        # Create GIF
        with imageio.get_writer('gif/response_matrix.gif', mode='I', duration=0.1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)
        for filename in filenames:
            os.remove(filename)

    def plotSingleResponse(R, true_energies):
        _, ax = plt.subplots(figsize=(8, 5))

        indices = [30, len(true_energies)//2, len(true_energies)-30]  # first, middle, last
        labels = ['Low edge', 'Mid range', 'High edge']

        for idx, label in zip(indices, labels):
            ax.plot(true_energies, R[:, idx], label=f"{label} ({true_energies[idx]:.2f} ph)")

        ax.set_xlabel("Photon Count")
        ax.set_ylabel("Probability Density")
        ax.set_title("Response Matrix Edge Cases")
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        plt.show()
    
    def saveResponseMatrix(R, output_filename="response_matrix.pdf", show=False):
        im = plt.imshow(R, aspect='auto', origin='lower',
                    extent=[true_energies[0], true_energies[-1],
                            measured_energies[0], measured_energies[-1]],
                    norm=LogNorm(vmin=1e-4, vmax=1))

        plt.xlabel('True Photon Count')
        plt.ylabel('Measured Photon Count')
        plt.title('Detector Response Matrix')
        plt.colorbar(im)
        plt.savefig(output_filename)
        if(show):
            plt.show()
        plt.close()
    
    def saveClosurePlot(unfolder, output_filename="chi2_convergence.pdf", show=False):
        plt.rcParams.update({'font.size': 32})
        plt.figure(figsize=(16,12))
        plt.plot(range(len(unfolder.chi2))[1:], unfolder.chi2[1:], 'bo-')
        plt.xlabel('Iteration')
        plt.ylabel('Chi-squared/N_DOF')
        plt.title('Convergence of Unfolding')
        plt.grid(True)

        plt.tight_layout()
        plt.savefig(output_filename, dpi=300)
        if(show):
            plt.show()
        plt.close()
    
    def saveUnfoldedComparison(true_energies, measured_spectrum, unfolded_spectrum, output_filename="unfoldedSpectrum.pdf", peaks=None, show=False):
        plt.rcParams.update({'font.size': 32})
        plt.figure(figsize=(16,12))
        plt.step(true_energies, measured_spectrum, color="#000000", linewidth=2, label='Measured Spectrum')
        plt.plot(true_energies, unfolded_spectrum, color=targetColorMap[0], linewidth=5, label=f'Unfolded Spectrum')
        linestyles = ["--", ":"]
        if(peaks):
            for idx, x_fit in enumerate(peaks[1]):
                param = peaks[0][idx]
                plt.plot(x_fit, gaussian(x_fit, *param[0]), targetColorMap[0], linewidth=10,linestyle=linestyles[idx], label = fr"Peak at {param[0][1]:.2f} ± {param[1][2]:.2f}")
        plt.xlabel('True Photon Count')
        plt.ylabel('Counts')
        # plt.title('Measured vs Unfolded Spectrum')
        plt.legend()
        plt.grid(True)
        plt.xlim(0, 350)
        plt.yscale("log")
        plt.ylim(bottom = 0.1)
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        if(show):
            plt.show()
        plt.close()
    
    def savePredictedComparison(measured_energies, measured_spectrum, predicted_measured, output_filename="predictedMeasuredSpectrum.pdf", show=False):
        plt.rcParams.update({'font.size': 32})
        plt.figure(figsize=(16,12))
        plt.step(measured_energies, measured_spectrum, color="#000000", linewidth=2, label='Measured Spectrum')
        plt.plot(measured_energies[1:], predicted_measured[1:], color=targetColorMap[0], linewidth=5, label=f'Predicted from Unfolded')
        plt.xlabel('Measured Photon Count')
        plt.ylabel('Counts')
        # plt.title('Measured vs Predicted Spectrum')
        plt.legend()
        plt.grid(True)
        plt.xlim(0, 350)
        plt.yscale("log")
        plt.ylim(bottom = 0.1)
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        if(show):
            plt.show()
        plt.close()
        
    print("Loading measured spectrum...")   

    mean_SPE = 110.342
    sigma_SPE= 10.0004
    mean_ped= 85.5295
    sigma_ped= 0.48595
    PDE = 0.2316
    saveMatrix = False
    debug = True
    lateral = True

    datapath = "../lightyield/data/2504/flat"
    if lateral:
        datapath = "../lightyield/data/2504/lateral"
    input_files = []
    output_files = []
    for crystal in range(36):
        output_files.append(f"{datapath}/output/crystal{crystal}")
        if lateral:
            input_files.append(f"{datapath}/bps_pwo_lateral_{crystal}_na22.his.txt")
        else:
            input_files.append(f"{datapath}/bps_pwo_{crystal}_na22.his.txt")
        
    measured_spectrum = []
    pe_axis = []

    with open(input_files[0], "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            x, y = map(float, parts[:2])
            pe_axis.append(x)
            measured_spectrum.append(y)

    pe_axis = np.array(pe_axis)
    pedastal_idx = np.argmin(np.abs(pe_axis - mean_ped))
    pe_axis = (pe_axis-mean_ped) / (mean_SPE-mean_ped)
    
    measured_spectrum = measured_spectrum[pedastal_idx:]
    measured_spectrum = np.array(measured_spectrum)


    unfolder = BayesianUnfoldingWithDRF(n_true_bins=4001-pedastal_idx, n_measured_bins=4001-pedastal_idx, n_iterations=40)

    maxE = (4001-pedastal_idx-mean_ped)/(mean_SPE-mean_ped)/PDE
    
    useExistingResponseMatrix = True
    if(not useExistingResponseMatrix):
        print("Building response matrix from DRF function...")
        R, true_energies, measured_energies = unfolder.build_response_matrix(
            true_energy_range=(0, maxE),      # Adjust to your energy range
            measured_energy_range=(0, maxE),
            custom_drf=pmt_response            # Use your DRF function
        )
    else:
        print("Loading response matrix from DRF function...")
        R, true_energies, measured_energies = unfolder.load_response_matrix(
            true_energy_range=(0, maxE), 
            measured_energy_range=(0, maxE),
            filename="responseMatrix.root")
        
    noisePeak = gaussian(true_energies, amplitude=3775.71, mean=17.9493, sigma=0.990767)
    # measured_spectrum = np.maximum(0, measured_spectrum - noisePeak)
    
    if (saveMatrix):
        responseShape  = R.shape
        hist1d = ROOT.TH1D("h_measurement", "Measurement", responseShape[0], 0, maxE)
        hist2d = ROOT.TH2D("h_responseMatrix", "Response Matrix", responseShape[0], 0, maxE, responseShape[1], 0, maxE)

        for i in range(responseShape[0]):
            for j in range(responseShape[1]):
                hist2d.SetBinContent(i+1, j+1, R[i, j])  # ROOT bins start at 1
        for i in range(responseShape[0]):
            hist1d.SetBinContent(i+1, measured_spectrum[i])  # ROOT bins start at 1

        root_file = ROOT.TFile("responseMatrix.root", "RECREATE")
        hist2d.Write()
        hist1d.Write()
        root_file.Close()

    print("\nPerforming Bayesian unfolding...")

    initial_prior = np.ones_like(true_energies)
    useold = True
    unfolded_spectrum = []
    predicted_measured = []
    if(useold):
        unfolded_spectrum = unfolder.unfold(measured_spectrum, prior=initial_prior)
        #unfolded_spectrum = unfolder.unfold_regularized(measured_spectrum, prior=measured_spectrum)
        predicted_measured = np.dot(unfolder.R, unfolded_spectrum)
    else:
        ch = 0
        initial_prior = np.ones_like(true_energies)
        true_edges = np.array(calcBinEdges(true_energies))
        measuredEdges = np.array(calcBinEdges(measured_energies))
        
        h_meas  = numpy_to_th1(f"h_meas{ch}", "Measured", measured_spectrum, measuredEdges)

        h_true  = numpy_to_th1(f"h_true{ch}", "Truth", initial_prior, true_edges)
        
        h_resp  = numpy_to_th2(f"h_resp{ch}", "Response", R, measuredEdges, true_edges)
        response = ROOT.RooUnfoldResponse(h_meas, h_true, h_resp)
        unfold = ROOT.RooUnfoldBayes(response , h_meas, 10)
        h_unfold = unfold.Hunfold()
        print(f"Unfolding complete of channel {ch}.")
        cov_matrix = unfold.Eunfold()
        print("Covariance matrix retrieved.")
        
        x_meas, y_meas = th1_to_numpy(h_meas)
        
        x_true, y_true = th1_to_numpy(h_true)
        x_unf, y_unf = th1_to_numpy(h_unfold)
            
        scale = np.sum(y_meas) / np.sum(y_unf) if np.sum(y_unf) !=0 else 1
        if(scale==0):
            scale=1
        print(f"Scale: {scale}")
        h_unfold.Scale(scale)
        x_unf, y_unf = th1_to_numpy(h_unfold)
        n_bins = h_unfold.GetNbinsX()
        cov = np.zeros((n_bins, n_bins))

        for i in range(n_bins):
            for j in range(n_bins):
                cov[i, j] = cov_matrix[i][j]*scale**2
        bin_errs = np.sqrt(np.diag(cov))


        plt.figure(figsize=(10,6))
        plt.step(measured_energies, measured_spectrum, where='mid', label='Measured', color='tab:orange', linewidth=1)
        plt.plot(x_unf, y_unf, label='Unfolded', color='tab:blue', linewidth=1)
        plt.fill_between(x_unf, y_unf - bin_errs, y_unf + bin_errs, color='tab:red', alpha=0.3)
        plt.xlabel('Deposited Energy / MeV')
        plt.ylabel('Counts')
        plt.title(f'Channel {ch}: Unfolded Spectrum vs Measured Energy Spectrum')
        plt.legend()
        plt.grid(True)
        plt.show()

        h_predicted = response.ApplyToTruth(h_unfold)# or .Hmeas()
        x_pred, y_pred = th1_to_numpy(h_predicted)
        
        plt.plot(x_pred, y_pred, '--', label='Predicted Measured', color='tab:blue', linewidth=1)
        plt.step(x_meas, y_meas, where='mid', label='Measured', color='tab:orange', linewidth=1)
        plt.xlabel('Measured ADC Counts')
        plt.ylabel('Counts')
        plt.title(f'Channel {ch}: Measured vs Predicted')
        plt.legend()
        plt.grid(True)
        plt.show()

    print("Plotting results...")
    if(saveMatrix):
        file = ROOT.TFile("my_hist.root", "RECREATE")

        hist = ROOT.TH1D("h_unfold", "Energy / PE", len(true_energies), true_energies[0], true_energies[-1])
        hist2 = ROOT.TH1D("h_measurement", "Energy / PE", len(true_energies), true_energies[0], true_energies[-1])
        hist3 = ROOT.TH1D("h_predicted", "Energy / PE", len(true_energies), true_energies[0], true_energies[-1])
        for index in range(len(true_energies)):
            hist.SetBinContent(index, unfolded_spectrum[index])
            hist2.SetBinContent(index, measured_spectrum[index])
            hist3.SetBinContent(index, predicted_measured[index])

        hist.Write()
        hist2.Write()
        hist3.Write()
        file.Close()
        
    selected_peaks, params, x_fits = findPeaks(unfolded_spectrum, measured_spectrum, true_energies)

    R[R < 1e-10] = 1e-10
    targetColorMap = ["#1f77b4", "#4e79a7", "#76b7b2", "#bab0ac", "#f28e2b", "#e15759", "#9c755f"]
    if (debug):
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        im = axes[0, 0].imshow(R, aspect='auto', origin='lower', 
                              extent=[true_energies[0], true_energies[-1], 
                                      measured_energies[0], measured_energies[-1]],
                                      norm=LogNorm(vmin=1e-4, vmax=1))
        axes[0, 0].set_xlabel('True Photon Count')
        axes[0, 0].set_ylabel('Measured Photon Count')
        axes[0, 0].set_title('Detector Response Matrix from DRF')
        plt.colorbar(im, ax=axes[0, 0], label='Probability')

        min_index = np.argmin(unfolder.chi2[1:]) + 1 


        axes[0, 1].step(true_energies, measured_spectrum, color="#000000", linewidth=2, label='Measured Spectrum')
        axes[0, 1].step(true_energies, unfolded_spectrum, color=targetColorMap[0], linewidth=2, label=f'Unfolded Spectrum')
        # axes[0, 1].plot(true_energies[peaks], unfolded_spectrum[peaks], 'rx', linewidth=2, label=f'Peaks of Unfolded Spectrum')
        for i, param in enumerate(params):
            x_fit = x_fits[i]
            # axes[0, 1].plot(x_fit, gaussian(x_fit, *param[0]), 'r--' , label = fr"Peak at {param[0][1]:.2f} ± {param[1][2]:.2f}")
        axes[0, 1].set_xlabel('Photon Count')
        axes[0, 1].set_ylabel('Counts')
        #axes[0, 1].set_yscale('log')
        axes[0, 1].set_title('True vs Unfolded Spectrum')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)

        axes[1, 0].step(measured_energies, measured_spectrum, color="#000000", linewidth=2, label='Measured Spectrum')
        axes[1, 0].step(measured_energies[1:], predicted_measured[1:], color=targetColorMap[0], linewidth=2, label=f'Predicted from Unfolded')
        # axes[1, 0].step(measured_energies, np.dot(unfolder.R, unfolder.history[min_index]), 'g-', linewidth=2, label=f'Unfolded Spectrum {min_index+1}')
        axes[1, 0].set_xlabel('Photon Count')
        axes[1, 0].set_ylabel('Counts')
        # axes[1, 0].set_yscale('log')
        axes[1, 0].set_title('Measured Spectrum vs Prediction')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)

        iterations = range(len(unfolder.chi2))
        axes[1, 1].plot(iterations, unfolder.chi2, 'bo-')
        # axes[1, 1].plot(iterations, unfolder.closure_errors, 'go-')
        axes[1, 1].set_xlabel('Iteration')
        axes[1, 1].set_ylabel('Chi-squared/N_DOF')
        axes[1, 1].set_title('Convergence of Unfolding')
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.close()
        # plt.show()

    #analyze all
    saveResponseMatrix(R, output_filename=f"{datapath}/output/responseMatrix.pdf")
    lightyields = []
    for crystalidx, infile in enumerate(input_files):
        print(f"\nAnalyzing crystal {crystalidx}...")
        measured_spectrum = []
        with open(infile, "r") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                _, y = map(float, parts[:2])
                measured_spectrum.append(y)
        measured_spectrum = measured_spectrum[pedastal_idx:]
        measured_spectrum = np.array(measured_spectrum)
        
        #initial_prior = np.ones_like(true_energies)
        ##unfold
        unfolded_spectrum = unfolder.unfold(measured_spectrum, prior=initial_prior)
        print(f"Chi^2: {unfolder.chi2[-1]}")
        predicted_measured = np.dot(unfolder.R, unfolded_spectrum)
        selected_peaks, params, x_fits = findPeaks(unfolded_spectrum, measured_spectrum, true_energies)
        
        print(f"selected_peaks {selected_peaks}")
        print(f"params {params}")
        print(f"measured_spectrum {measured_spectrum[selected_peaks]}")
        print(f"true_energies {true_energies[selected_peaks]}")
        
        y1 = params[0][0][1]
        y2 = params[1][0][1]
        
        x1 = 0.51099
        x2 = 1.27453
        
        sy1 = params[0][1][1]
        sy2 = params[1][1][1]
        
        print(f"x: {x1}; {x2}")
        print(f"y: {y1}; {y2}")
        print(f"sy: {sy1}; {sy2}")
        
        m = (y2 - y1) / (x2 - x1)

        sm = np.sqrt(sy1**2 + sy2**2) / abs(x2 - x1)**2
        
        print(f"m: {m}; sm {sm}")
        lightyields.append([m, sm])
        saveUnfoldedComparison(true_energies, measured_spectrum, unfolded_spectrum, output_filename=f"{output_files[crystalidx]}/unfoldedSpectrum.pdf", peaks=[params, x_fits], show=False)
        savePredictedComparison(measured_energies, measured_spectrum, predicted_measured, output_filename=f"{output_files[crystalidx]}/predictedSpectrum.pdf", show=False)
        saveClosurePlot(unfolder, output_filename=f"{output_files[crystalidx]}/chi2_Convergence.pdf", show=False)
    
    np.save(f'{datapath}/output/fit_params.npy', np.array(lightyields, dtype=object))

    

