import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import curve_fit

from scipy.stats import poisson
from scipy.signal import convolve
import matplotlib.gridspec as gridspec

class BayesianUnfoldingWithDRF:
    """
    Bayesian Unfolding with a functional Detector Response Function (DRF)
    """
    
    def __init__(self, n_true_bins=50, n_measured_bins=50, n_iterations=10):
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
        
        # Build response matrix
        response_matrix = np.zeros((self.n_measured_bins, self.n_true_bins))
        
        print("Building response matrix...")
        for j, true_energy in enumerate(self.true_centers):
            if j % 10 == 0:
                print(f"  Processing true energy bin {j+1}/{self.n_true_bins}")
            
            # For each true energy, calculate probability in each measured bin
            for i in range(self.n_measured_bins):
                # Integrate DRF over the measured energy bin
                low = self.measured_edges[i]
                high = self.measured_edges[i + 1]
                
                # Numerical integration of DRF over the bin
                integral, error = quad(drf, low, high, args=(true_energy,))
                response_matrix[i, j] = integral
        
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
            f_new = self._smooth_spectrum(f_new, strength=0.1)
            
            # Ensure non-negativity
            f_new = np.maximum(f_new, 0)
            
            f = f_new
            self.history.append(f.copy())
            
            print(f"Iteration {iteration+1}: Total counts = {np.sum(f):.1f}")
        
        self.unfolded = f
        return f
    
    def _smooth_spectrum(self, spectrum, strength=0.1):
        """Apply mild smoothing to reduce oscillations"""
        from scipy.ndimage import gaussian_filter1d
        return (1 - strength) * spectrum + strength * gaussian_filter1d(spectrum, sigma=1.0)
    
    def calculate_closure(self, measured_spectrum):
        """Calculate how well unfolded spectrum reproduces measurement"""
        if hasattr(self, 'unfolded'):
            predicted_measured = np.dot(self.R, self.unfolded)
            closure = np.sum((measured_spectrum - predicted_measured)**2) / np.sum(measured_measured**2)
            return closure, predicted_measured
        else:
            raise ValueError("Must run unfold() first")

# Example usage with your data
if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(42)
    
    print("Bayesian Unfolding with Functional DRF")
    print("=" * 50)
    
    # 1. Initialize unfolding with your binning preferences
    unfolder = BayesianUnfoldingWithDRF(n_true_bins=50, n_measured_bins=50, n_iterations=15)
    
    # 2. DEFINE YOUR ACTUAL DRF FUNCTION HERE
    def my_actual_drf(measured_energy, true_energy):
        """
        REPLACE THIS WITH YOUR ACTUAL DETECTOR RESPONSE FUNCTION!
        
        Parameters:
        measured_energy: the energy you measure (in your units)
        true_energy: the true deposited energy (in your units)
        
        Returns:
        probability density for measuring 'measured_energy' given 'true_energy'
        """
        # Example: Gaussian response with energy-dependent resolution
        # Modify this to match your actual detector response
        
        # Your energy resolution (from your light yield measurements)
        # This is just an example - use your actual resolution function
        light_yield = 5000  # photons/MeV (replace with your value)
        n_pe = true_energy * light_yield  # mean number of photoelectrons
        sigma_pe = np.sqrt(n_pe)  # standard deviation from Poisson statistics
        
        # Convert to energy resolution
        sigma_energy = sigma_pe / light_yield
        
        # Gaussian response
        return np.exp(-0.5 * ((measured_energy - true_energy) / sigma_energy)**2) / (sigma_energy * np.sqrt(2 * np.pi))
    
    # 3. Build response matrix using your DRF
    print("Building response matrix from DRF function...")
    R, true_energies, measured_energies = unfolder.build_response_matrix(
        true_energy_range=(0.1, 5.0),      # Adjust to your energy range
        measured_energy_range=(0.1, 5.0),   # Adjust to your energy range  
        custom_drf=my_actual_drf            # Use your DRF function
    )
    
    # 4. CREATE OR LOAD YOUR MEASURED SPECTRUM
    # Replace this with your actual measured data
    print("Loading measured spectrum...")
    
    # Example: create a test measured spectrum
    def test_true_spectrum(energy):
        """Create a test true spectrum for demonstration"""
        return (np.exp(-energy/2.0) + 
                0.5 * np.exp(-0.5 * ((energy - 1.5) / 0.1)**2) + 
                0.3 * np.exp(-0.5 * ((energy - 2.5) / 0.08)**2))
    
    true_spec = test_true_spectrum(true_energies) * 10000
    expected_measured = np.dot(R, true_spec)
    measured_spectrum = np.random.poisson(expected_measured)  # Add Poisson noise
    
    # If you have your own measured data, replace the line above with:
    # measured_spectrum = your_actual_measured_data_array
    
    # 5. Perform unfolding
    print("\nPerforming Bayesian unfolding...")
    
    # You can provide a prior if you have knowledge about the expected spectrum
    # If not, a flat prior will be used
    initial_prior = np.ones_like(true_energies)  # Flat prior
    
    unfolded_spectrum = unfolder.unfold(measured_spectrum, prior=initial_prior)
    
    # 6. Calculate closure test
    closure_error, predicted_measured = unfolder.calculate_closure(measured_spectrum)
    print(f"\nClosure test error: {closure_error:.6f}")
    
    # 7. Plot results
    print("Plotting results...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Response matrix
    im = axes[0, 0].imshow(R, aspect='auto', origin='lower', 
                          extent=[true_energies[0], true_energies[-1], 
                                  measured_energies[0], measured_energies[-1]])
    axes[0, 0].set_xlabel('True Energy (MeV)')
    axes[0, 0].set_ylabel('Measured Energy (MeV)')
    axes[0, 0].set_title('Detector Response Matrix from DRF')
    plt.colorbar(im, ax=axes[0, 0], label='Probability')
    
    # Plot 2: True vs Unfolded spectrum
    axes[0, 1].plot(true_energies, true_spec, 'b-', linewidth=2, label='True Spectrum')
    axes[0, 1].plot(true_energies, unfolded_spectrum, 'r--', linewidth=2, label='Unfolded Spectrum')
    axes[0, 1].set_xlabel('Energy (MeV)')
    axes[0, 1].set_ylabel('Counts')
    axes[0, 1].set_title('True vs Unfolded Spectrum')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Measured spectrum and prediction from unfolded
    axes[1, 0].plot(measured_energies, measured_spectrum, 'ko', markersize=3, label='Measured')
    axes[1, 0].plot(measured_energies, predicted_measured, 'r-', linewidth=2, label='Predicted from Unfolded')
    axes[1, 0].set_xlabel('Measured Energy (MeV)')
    axes[1, 0].set_ylabel('Counts')
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