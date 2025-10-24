import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import curve_fit

from scipy.stats import poisson
from scipy.signal import convolve
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import ROOT

class BayesianUnfoldingWithDRF:
    """
    Bayesian Unfolding with a functional Detector Response Function (DRF)
    """
    
    def __init__(self, n_true_bins=4001, n_measured_bins=4001, n_iterations=10):
        self.n_true_bins = n_true_bins
        self.n_measured_bins = n_measured_bins
        self.n_iterations = n_iterations
    
    def pmt_response(true_adc, mean_SPE, sigma_SPE, mean_ped, sigma_ped, pe_range=None, n_max=300):
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

        mu_pe = (true_adc-mean_ped) / (mean_SPE-mean_ped)/0.2316

        if pe_range is None:
            pe_min = 0
            pe_max = (mu_pe + 7 * sigma_SPE)
            n_points = 1000
            pe_axis = np.linspace(pe_min, pe_max, n_points)
        else:
            pe_axis = np.linspace(pe_range[0], pe_range[1], pe_range[2])

        response = np.zeros_like(pe_axis)
        components = []

        for n in range(0, n_max + 1):    
            poisson_prob = poisson.pmf(n, mu_pe)
            if poisson_prob > 1e-10:  
                mean_n = n
                sigma_n = np.sqrt((np.sqrt(n) * sigma_SPE)**2+sigma_ped**2) if n > 0 else np.sqrt(sigma_SPE**2+sigma_ped**2)

                if sigma_n < 1e-10:
                    sigma_n = sigma_SPE

                gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                          np.exp(-0.5 * ((pe_axis - mean_n) / sigma_n) ** 2)

                component = poisson_prob * gaussian
                response += component
                components.append((n, mean_n, sigma_n, component))

        return pe_axis, response, components

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
        
        self.history = [f.copy()]
        
        print("\nStarting Bayesian unfolding...")
        for iteration in range(self.n_iterations):
            g_expected = np.dot(self.R, f)
            
            g_expected[g_expected == 0] = 1e-10
            
            correction_factor = np.dot(self.R.T, g / g_expected)
            f_new = f * correction_factor / efficiency
            
            f_new = self._smooth_spectrum(f_new, strength=0.1)
            
            f_new = np.maximum(f_new, 0)
            
            f = f_new
            self.history.append(f.copy())
            
            print(f"Iteration {iteration+1}: Total counts = {np.sum(f):.1f}")
        
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
        return (1 - strength) * spectrum + strength * gaussian_filter1d(spectrum, sigma=1.0)
    
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

# Example usage with your data
if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(42)
    
    print("Bayesian Unfolding with Functional DRF")
    print("=" * 50)

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

        for n in range(0, n_max + 1):    
            poisson_prob = poisson.pmf(n, true_pe*PDE)
            if poisson_prob > 1e-10:  
                sigma_n = np.sqrt((np.sqrt(n) * sigma_SPE)**2+sigma_ped**2) if n > 0 else np.sqrt(sigma_SPE**2+sigma_ped**2)

                if sigma_n < 1e-10:
                    sigma_n = sigma_SPE

                gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                          np.exp(-0.5 * ((pe_axis*PDE - n) / sigma_n) ** 2)

                component = poisson_prob * gaussian
                response += component

        #return pe_axis, response, components
        return response

    def pmt_response_single(measured_pe, true_pe):
        """
        Calculate PMT response function for a given energy deposition

        Parameters:
        - measured_pe: Measured pe value
        - true_pe: True pe value

        Returns:
        - response: Probability density function
        - components: Individual photoelectron peaks for plotting
        """
        n_max=50
        mean_SPE = 110.342
        sigma_SPE= 10.0004
        mean_ped= 85.5295
        sigma_ped= 0.48595
        sigma_SPE = (sigma_SPE) / (mean_SPE-mean_ped)
        sigma_ped = (sigma_ped) / (mean_SPE-mean_ped)
        components = []
        response = 0

        for n in range(0, n_max + 1):    
            poisson_prob = poisson.pmf(n, true_pe)
            if poisson_prob > 1e-10:  
                mean_n = n
                sigma_n = np.sqrt((np.sqrt(n) * sigma_SPE)**2+sigma_ped**2) if n > 0 else np.sqrt(sigma_SPE**2+sigma_ped**2)

                if sigma_n < 1e-10:
                    sigma_n = sigma_SPE

                gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                          np.exp(-0.5 * ((measured_pe - mean_n) / sigma_n) ** 2)

                component = poisson_prob * gaussian
                response += component
                components.append((n, mean_n, sigma_n, component))

        return response

    print("Loading measured spectrum...")
    
    
    # datapath = "../lightyield/data/2504/flat"
    # input_file = f"{datapath}/bps_pwo_3_na22.his.txt"
    
    datapath = "../lightyield/old"
    input_file = f"{datapath}/BPS_CH0.his.txt"
    
    mean_SPE = 110.342
    sigma_SPE= 10.0004
    mean_ped= 85.5295
    sigma_ped= 0.48595
    PDE = 0.2316

    print(f"Mean SPE: {1}, Sigma SPE: {sigma_SPE/mean_SPE}")
    print(f"Mean PED: {mean_ped/mean_SPE}, Sigma PED: {sigma_ped/mean_SPE}")
    # --- Read data from text file ---
    measured_spectrum = []
    pe_axis = []
    with open(input_file, "r") as f:
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
    unfolder = BayesianUnfoldingWithDRF(n_true_bins=4001-pedastal_idx, n_measured_bins=4001-pedastal_idx, n_iterations=10)
    # 3. Build response matrix using your DRF
    print("Building response matrix from DRF function...")

    maxE = (4001-pedastal_idx-mean_ped)/(mean_SPE-mean_ped)/PDE
    R, true_energies, measured_energies = unfolder.build_response_matrix(
        true_energy_range=(0, maxE),      # Adjust to your energy range
        measured_energy_range=(0, maxE),
        custom_drf=pmt_response            # Use your DRF function
    )
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

    initial_prior = np.ones_like(true_energies)  # Flat prior
    
    # unfolder.unfold_bayesian(measured_spectrum, measured_spectrum)
    unfolded_spectrum = unfolder.unfold(measured_spectrum, prior=measured_spectrum)
    #unfolded_spectrum = unfolder.unfold_regularized(measured_spectrum, prior=measured_spectrum)
    
    # 6. Calculate closure test
    closure_error, chi2_per_dof, predicted_measured = unfolder.calculate_closure(measured_spectrum)
    print(f"Closure test error: {closure_error:.4f} ({closure_error*100:.1f}%)")
    print(f"Chi-squared per DOF: {chi2_per_dof:.2f}")
    
    closure_results = unfolder.calculate_closure_fixed(measured_spectrum)
    predicted_measured = closure_results['predicted']

    print(f"\n=== INTERPRETATION GUIDE ===")
    print(f"RMS Error < 10%: Excellent")
    print(f"RMS Error 10-20%: Good") 
    print(f"RMS Error 20-30%: Acceptable")
    print(f"RMS Error > 30%: Needs improvement")
    print(f"Correlation > 0.95: Excellent")
    print(f"Correlation 0.9-0.95: Good")
    print(f"Correlation < 0.9: Needs improvement")

    unfolder.debug_closure_issue(measured_spectrum, predicted_measured)

    # 7. Plot results
    print("Plotting results...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    

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

    # Plot 1: Response matrix
    im = axes[0, 0].imshow(R, aspect='auto', origin='lower', 
                          extent=[true_energies[0], true_energies[-1], 
                                  measured_energies[0], measured_energies[-1]],
                                  norm=LogNorm(vmin=1e-4, vmax=1))
    axes[0, 0].set_xlabel('True Photoelectrons (p.e.)')
    axes[0, 0].set_ylabel('Measured Photoelectrons (p.e.)')
    axes[0, 0].set_title('Detector Response Matrix from DRF')
    plt.colorbar(im, ax=axes[0, 0], label='Probability')
    
    # Plot 2: True vs Unfolded spectrum
    axes[0, 1].step(true_energies, measured_spectrum, 'b-', linewidth=2, label='Measured Spectrum')
    axes[0, 1].step(true_energies, unfolded_spectrum, 'r--', linewidth=2, label='Unfolded Spectrum')
    axes[0, 1].set_xlabel('Photoelectrons (p.e.)')
    axes[0, 1].set_ylabel('Counts')
    #axes[0, 1].set_yscale('log')
    axes[0, 1].set_title('True vs Unfolded Spectrum')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Measured spectrum and prediction from unfolded
    axes[1, 0].step(measured_energies, measured_spectrum, 'ko', markersize=3, label='Measured')
    axes[1, 0].step(measured_energies, predicted_measured, 'r-', linewidth=2, label='Predicted from Unfolded')
    axes[1, 0].set_xlabel('Measured Photoelectrons (p.e.)')
    axes[1, 0].set_ylabel('Counts')
    # axes[1, 0].set_yscale('log')
    axes[1, 0].set_title('Measured Spectrum vs Prediction')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Convergence
    iterations = range(len(unfolder.history))
    total_counts = [np.sum(f) for f in unfolder.history]
    axes[1, 1].plot(iterations, total_counts, 'bo-')
    axes[1, 1].set_xlabel('Iteration')
    axes[1, 1].set_ylabel('Total Counts in Unfolded Spectrum')
    axes[1, 1].set_title('Convergence of Unfolding')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    print("\nUnfolding completed successfully!")