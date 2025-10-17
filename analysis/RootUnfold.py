import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import ROOT
from array import array

from scipy.integrate import quad
from scipy.optimize import curve_fit

from scipy.stats import poisson
from scipy.signal import convolve
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm


class RooUnfoldCrossCheck:
    """
    Cross-check our Bayesian unfolding with RooUnfold
    """
    
    def __init__(self):
        # Make sure ROOT and RooUnfold are available
        try:
            ROOT.gSystem.Load("libRooUnfold")
            print("RooUnfold loaded successfully")
        except:
            print("Warning: RooUnfold not available. Some features will be disabled.")
    
    def convert_to_root_histograms(self, true_edges, measured_edges, response_matrix, 
                                 measured_spectrum, true_spectrum=None):
        """
        Convert numpy arrays to ROOT histograms for RooUnfold
        """
        # Convert edges to arrays for ROOT
        true_edges_root = array('d', true_edges)
        measured_edges_root = array('d', measured_edges)
        
        n_true_bins = len(true_edges) - 1
        n_measured_bins = len(measured_edges) - 1
        
        # Create histograms
        h_measured = ROOT.TH1D("h_measured", "Measured Spectrum", 
                              n_measured_bins, measured_edges_root)
        h_true = None
        if true_spectrum is not None:
            h_true = ROOT.TH1D("h_true", "True Spectrum", 
                              n_true_bins, true_edges_root)
        
        # Create response matrix in ROOT format
        h_response = ROOT.TH2D("h_response", "Response Matrix",
                              n_measured_bins, measured_edges_root,
                              n_true_bins, true_edges_root)
        
        # Fill measured spectrum
        for i in range(n_measured_bins):
            h_measured.SetBinContent(i+1, measured_spectrum[i])
            h_measured.SetBinError(i+1, np.sqrt(measured_spectrum[i] + 1e-10))
        
        # Fill true spectrum if available
        if true_spectrum is not None:
            for j in range(n_true_bins):
                h_true.SetBinContent(j+1, true_spectrum[j])
        
        # Fill response matrix
        for i in range(n_measured_bins):
            for j in range(n_true_bins):
                h_response.SetBinContent(i+1, j+1, response_matrix[i, j])
        
        return h_measured, h_true, h_response
    
    def run_roounfold_bayes(self, true_edges, measured_edges, response_matrix,
                           measured_spectrum, true_spectrum=None, n_iterations=4):
        """
        Run Bayesian unfolding using RooUnfold
        """
        print("\n=== RUNNING ROOUNFOLD BAYESIAN UNFOLDING ===")
        
        # Convert to ROOT histograms
        h_meas, h_true, h_resp = self.convert_to_root_histograms(
            true_edges, measured_edges, response_matrix, measured_spectrum, true_spectrum)
        
        # Create RooUnfold response object
        response = ROOT.RooUnfoldResponse(h_meas, h_true, h_resp)
        
        # Create and run Bayesian unfolding
        unfold_bayes = ROOT.RooUnfoldBayes(response, h_meas, n_iterations)
        h_unfolded = unfold_bayes.Hreco()
        
        # Get uncertainty
        h_errors = unfold_bayes.Ereco()
        
        # Convert back to numpy
        n_true_bins = len(true_edges) - 1
        unfolded_roounfold = np.zeros(n_true_bins)
        errors_roounfold = np.zeros(n_true_bins)
        
        for i in range(n_true_bins):
            unfolded_roounfold[i] = h_unfolded.GetBinContent(i+1)
            errors_roounfold[i] = h_errors.GetBinContent(i+1)
        
        # Get chi-squared and other diagnostics
        chi2 = unfold_bayes.Chi2()
        dof = n_true_bins
        
        print(f"RooUnfold Bayes completed:")
        print(f"  Iterations: {n_iterations}")
        print(f"  Chi-squared: {chi2:.2f}")
        print(f"  Total unfolded counts: {np.sum(unfolded_roounfold):.1f}")
        
        return unfolded_roounfold, errors_roounfold, chi2
    
    def run_roounfold_svd(self, true_edges, measured_edges, response_matrix,
                         measured_spectrum, true_spectrum=None, k_reg=3):
        """
        Run SVD unfolding using RooUnfold (alternative method)
        """
        print("\n=== RUNNING ROOUNFOLD SVD UNFOLDING ===")
        
        h_meas, h_true, h_resp = self.convert_to_root_histograms(
            true_edges, measured_edges, response_matrix, measured_spectrum, true_spectrum)
        
        response = ROOT.RooUnfoldResponse(h_meas, h_true, h_resp)
        
        # SVD unfolding
        unfold_svd = ROOT.RooUnfoldSvd(response, h_meas, k_reg)
        h_unfolded = unfold_svd.Hreco()
        h_errors = unfold_svd.Ereco()
        
        n_true_bins = len(true_edges) - 1
        unfolded_svd = np.zeros(n_true_bins)
        errors_svd = np.zeros(n_true_bins)
        
        for i in range(n_true_bins):
            unfolded_svd[i] = h_unfolded.GetBinContent(i+1)
            errors_svd[i] = h_errors.GetBinContent(i+1)
        
        print(f"RooUnfold SVD completed:")
        print(f"  Regularization parameter: {k_reg}")
        print(f"  Total unfolded counts: {np.sum(unfolded_svd):.1f}")
        
        return unfolded_svd, errors_svd
    
    def compare_methods(self, true_edges, measured_edges, response_matrix,
                       measured_spectrum, true_spectrum, our_unfolded,
                       our_closure_error, n_iterations=4):
        """
        Compare our implementation with RooUnfold
        """
        print("\n" + "="*60)
        print("METHOD COMPARISON: Our Implementation vs RooUnfold")
        print("="*60)
        
        # Run RooUnfold Bayesian unfolding
        unfolded_roo, errors_roo, chi2_roo = self.run_roounfold_bayes(
            true_edges, measured_edges, response_matrix, 
            measured_spectrum, true_spectrum, n_iterations)
        
        # Calculate metrics for both methods
        true_centers = (true_edges[1:] + true_edges[:-1]) / 2
        
        # Calculate closure for RooUnfold
        predicted_measured_roo = np.dot(response_matrix, unfolded_roo)
        mask = measured_spectrum > np.max(measured_spectrum) * 0.01
        if np.sum(mask) > 0:
            closure_roo = np.sqrt(np.mean(
                ((measured_spectrum[mask] - predicted_measured_roo[mask]) / 
                 measured_spectrum[mask])**2))
        else:
            closure_roo = float('inf')
        
        # Calculate correlation between methods
        correlation = np.corrcoef(our_unfolded, unfolded_roo)[0, 1]
        
        # Calculate relative difference
        relative_diff = np.abs(our_unfolded - unfolded_roo) / (unfolded_roo + 1e-10)
        mean_relative_diff = np.mean(relative_diff)
        
        print(f"\n=== COMPARISON RESULTS ===")
        print(f"Our Implementation:")
        print(f"  Total counts: {np.sum(our_unfolded):.1f}")
        print(f"  Closure error: {our_closure_error:.4f} ({our_closure_error*100:.1f}%)")
        print(f"RooUnfold:")
        print(f"  Total counts: {np.sum(unfolded_roo):.1f}")
        print(f"  Closure error: {closure_roo:.4f} ({closure_roo*100:.1f}%)")
        print(f"  Chi-squared: {chi2_roo:.2f}")
        print(f"Comparison:")
        print(f"  Correlation between methods: {correlation:.4f}")
        print(f"  Mean relative difference: {mean_relative_diff:.4f} ({mean_relative_diff*100:.1f}%)")
        
        # Plot comparison
        self.plot_comparison(true_centers, true_spectrum, our_unfolded, 
                           unfolded_roo, errors_roo, correlation, mean_relative_diff)
        
        return unfolded_roo, errors_roo, correlation, mean_relative_diff
    
    def plot_comparison(self, true_energies, true_spectrum, our_unfolded,
                       roounfold_unfolded, roounfold_errors, correlation, mean_diff):
        """Plot comparison between our method and RooUnfold"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Spectrum comparison
        axes[0, 0].plot(true_energies, true_spectrum, 'k-', linewidth=2, label='True Spectrum')
        axes[0, 0].plot(true_energies, our_unfolded, 'b--', linewidth=2, label='Our Unfolding')
        axes[0, 0].plot(true_energies, roounfold_unfolded, 'r:', linewidth=2, label='RooUnfold')
        axes[0, 0].fill_between(true_energies, 
                               roounfold_unfolded - roounfold_errors,
                               roounfold_unfolded + roounfold_errors,
                               alpha=0.3, color='red', label='RooUnfold Errors')
        axes[0, 0].set_xlabel('True Energy')
        axes[0, 0].set_ylabel('Counts')
        axes[0, 0].set_title('Unfolding Comparison')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Ratio to true spectrum
        ratio_ours = our_unfolded / (true_spectrum + 1e-10)
        ratio_roo = roounfold_unfolded / (true_spectrum + 1e-10)
        
        axes[0, 1].plot(true_energies, ratio_ours, 'b--', linewidth=2, label='Our Method / True')
        axes[0, 1].plot(true_energies, ratio_roo, 'r:', linewidth=2, label='RooUnfold / True')
        axes[0, 1].axhline(y=1.0, color='k', linestyle='-', alpha=0.5)
        axes[0, 1].set_xlabel('True Energy')
        axes[0, 1].set_ylabel('Ratio to True')
        axes[0, 1].set_title('Ratio to True Spectrum')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_ylim(0, 2)
        
        # Plot 3: Difference between methods
        relative_diff = (our_unfolded - roounfold_unfolded) / (roounfold_unfolded + 1e-10)
        axes[1, 0].plot(true_energies, relative_diff * 100, 'g-', linewidth=2)
        axes[1, 0].axhline(y=0, color='k', linestyle='-', alpha=0.5)
        axes[1, 0].set_xlabel('True Energy')
        axes[1, 0].set_ylabel('Relative Difference (%)')
        axes[1, 0].set_title(f'Our Method vs RooUnfold\nMean Difference: {mean_diff*100:.1f}%')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Correlation scatter plot
        axes[1, 1].scatter(our_unfolded, roounfold_unfolded, alpha=0.6)
        min_val = min(np.min(our_unfolded), np.min(roounfold_unfolded))
        max_val = max(np.max(our_unfolded), np.max(roounfold_unfolded))
        axes[1, 1].plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8)
        axes[1, 1].set_xlabel('Our Unfolded Counts')
        axes[1, 1].set_ylabel('RooUnfold Counts')
        axes[1, 1].set_title(f'Correlation: {correlation:.4f}')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()


class BayesianUnfoldingWithDRF:
    """
    Bayesian Unfolding with a functional Detector Response Function (DRF)
    """
    
    def __init__(self, n_true_bins=4001, n_measured_bins=4001, n_iterations=10):
        self.n_true_bins = n_true_bins
        self.n_measured_bins = n_measured_bins
        self.n_iterations = n_iterations
    
    def pmt_response(true_adc, mean_SPE, sigma_SPE, mean_ped, sigma_ped, pe_range=None, n_max=50):
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

        mu_pe = (true_adc-mean_ped) / (mean_SPE-mean_ped)

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
        
        # Estimate efficiency from response matrix if not provided
        if efficiency is None:
            efficiency = np.sum(self.R, axis=0)
        
        # Store iteration history
        self.history = [f.copy()]
        
        print("\nStarting Bayesian unfolding...")
        for iteration in range(self.n_iterations):
            # Calculate expected measured spectrum given current truth estimate
            g_expected = np.dot(self.R, f)
            
            # Avoid division by zero
            g_expected[g_expected == 0] = 1e-10
            
            # Bayesian update
            correction_factor = np.dot(self.R.T, g / g_expected)
            f_new = f * correction_factor / efficiency
            
            # Apply smoothing to reduce oscillations
            #f_new = self._smooth_spectrum(f_new, strength=0.1)
            
            # Ensure non-negativity
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
            # Expected measured spectrum
            g_expected = np.dot(self.R, f)
            g_expected[g_expected == 0] = 1e-10

            # Bayesian update
            correction_factor = np.dot(self.R.T, g / g_expected)
            f_new = f * correction_factor / efficiency

            # Apply Tikhonov regularization
            f_new = self._tikhonov_regularize(f_new, f, regularization_strength)

            # Ensure non-negativity
            f_new = np.maximum(f_new, 0)

            # Check convergence
            relative_change = np.sum(np.abs(f_new - f)) / np.sum(f)

            f = f_new
            self.history.append(f.copy())

            # Calculate closure error for monitoring
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

def pmt_response(true_pe, pe_axis, n_max=300):
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
    sigma_SPE = (sigma_SPE) / (mean_SPE-mean_ped)
    sigma_ped = (sigma_ped) / (mean_SPE-mean_ped)
    response = np.zeros_like(pe_axis)
    for n in range(0, n_max + 1):    
        poisson_prob = poisson.pmf(n, true_pe)
        if poisson_prob > 1e-10:  
            mean_n = n
            sigma_n = np.sqrt((np.sqrt(n) * sigma_SPE)**2+sigma_ped**2) if n > 0 else np.sqrt(sigma_SPE**2+sigma_ped**2)
            if sigma_n < 1e-10:
                sigma_n = sigma_SPE
            gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                      np.exp(-0.5 * ((pe_axis - mean_n) / sigma_n) ** 2)
            component = poisson_prob * gaussian
            response += component
    #return pe_axis, response, components
    return response

# Complete cross-check workflow
def full_cross_check():
    """Complete workflow for cross-checking our implementation"""
    
    # Initialize both implementations
    our_unfolder = BayesianUnfoldingWithDRF(n_true_bins=4001-85, n_measured_bins=4001-85, n_iterations=50)
    roo_checker = RooUnfoldCrossCheck()
    
    # Build response matrix (using your DRF)
    R, true_edges, measured_edges = our_unfolder.build_response_matrix(
        true_energy_range=(0, 157.76203526448361),
        measured_energy_range=(0, 157.76203526448361),
        custom_drf=pmt_response  # Your actual DRF
    )
    
    # Create or load your measured spectrum
    # measured_spectrum = your_actual_measured_data
    # true_spectrum = your_actual_true_spectrum (if available for testing)
    
    datapath = "../lightyield/data/2504/flat"
    input_file = f"{datapath}/bps_pwo_3_na22.his.txt"
    
    #datapath = "../lightyield/data/1310_pedastal"
    #input_file = f"{datapath}/bps_spe_100ns.his.txt"
    
    mean_SPE = 110.342
    sigma_SPE= 10.0004
    mean_ped= 85.5295
    sigma_ped= 0.48595

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
            pe_axis.append((x-mean_ped) / (mean_SPE-mean_ped))
            measured_spectrum.append(y)

    # For demonstration, create synthetic data
    true_centers = (true_edges[1:] + true_edges[:-1]) / 2
    def test_spectrum(e):
        return (np.exp(-e/2.0) + 
                0.5 * np.exp(-0.5 * ((e - 1.5) / 0.1)**2) +
                0.3 * np.exp(-0.5 * ((e - 2.5) / 0.08)**2))
    
    measured_spectrum = measured_spectrum[85:]
    print(len(measured_spectrum))
    measured_spectrum = np.array(measured_spectrum)
    pe_axis = np.array(pe_axis)
    
    # Run our implementation
    print("Running our Bayesian unfolding implementation...")
    our_unfolded = our_unfolder.unfold(measured_spectrum, prior=measured_spectrum)
    
    # Calculate closure for our method
    predicted_measured_ours = np.dot(R, our_unfolded)

    mask = measured_spectrum > np.max(measured_spectrum) * 0.01
    our_closure = np.sqrt(np.mean(
        ((measured_spectrum[mask] - predicted_measured_ours[mask]) / 
         measured_spectrum[mask])**2))
    
    # Cross-check with RooUnfold
    roo_unfolded, roo_errors, correlation, mean_diff = roo_checker.compare_methods(
        true_edges, measured_edges, R, measured_spectrum, measured_spectrum,
        our_unfolded, our_closure, n_iterations=10)
    
    # Print final assessment
    print("\n" + "="*60)
    print("FINAL ASSESSMENT")
    print("="*60)
    
    if correlation > 0.95:
        print("✅ EXCELLENT AGREEMENT: Our implementation matches RooUnfold very well!")
    elif correlation > 0.9:
        print("✅ GOOD AGREEMENT: Our implementation is consistent with RooUnfold")
    elif correlation > 0.8:
        print("⚠️  ACCEPTABLE AGREEMENT: Some differences, but generally consistent")
    else:
        print("❌ POOR AGREEMENT: Significant differences detected")
    
    if mean_diff < 0.1:
        print("✅ SMALL DIFFERENCES: Mean relative difference < 10%")
    elif mean_diff < 0.2:
        print("⚠️  MODERATE DIFFERENCES: Mean relative difference 10-20%")
    else:
        print("❌ LARGE DIFFERENCES: Mean relative difference > 20%")
    
    return our_unfolded, roo_unfolded, correlation, mean_diff

# Run the cross-check
if __name__ == "__main__":
    our_result, roo_result, corr, diff = full_cross_check()