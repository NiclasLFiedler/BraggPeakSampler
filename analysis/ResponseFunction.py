import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.signal import convolve
import matplotlib.gridspec as gridspec


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
    underflow_axis = np.arange(-10, 0, 0.17025946972945388)

    underflow_response = np.zeros_like(underflow_axis)
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
            
            if(true_pe<20):
                underflow_gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                          np.exp(-0.5 * ((underflow_axis*PDE - n) / sigma_n) ** 2)
                underflow_component = poisson_prob * underflow_gaussian

                underflow_response += underflow_component
                component[0] += np.sum(underflow_component)
            response += component
            components.append((n, n, sigma_n, component))
    # print()
    # print(response)
    return pe_axis, response, components
    #return response

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

def sipm_response(energy_deposited, light_yield, Q_pe, sigma_SPE, charge_range=None, n_max=50):
    """
    Calculate SiPM response function for a given energy deposition
    
    Parameters:
    - energy_deposited: True energy deposited (MeV)
    - light_yield: Light yield (photoelectrons/MeV)
    - Q_pe: Charge per photoelectron (arbitrary units)
    - sigma_SPE: Standard deviation of single photoelectron peak
    - charge_range: Tuple (min, max, n_points) for charge axis
    - n_max: Maximum number of photoelectrons to consider
    
    Returns:
    - charge_axis: Array of charge values
    - response: Probability density function
    - components: Individual photoelectron peaks for plotting
    """
    
    # Calculate mean number of photoelectrons
    mu_pe = energy_deposited * light_yield # * eta_PDE * eta_coll  (0.63, 0.2)
    
    # Set up charge axis if not provided
    if charge_range is None:
        charge_min = 0
        charge_max = (mu_pe + 4 * np.sqrt(mu_pe)) * Q_pe + 5 * sigma_SPE
        n_points = 1000
        charge_axis = np.linspace(charge_min, charge_max, n_points)
    else:
        charge_axis = np.linspace(charge_range[0], charge_range[1], charge_range[2])
    # Initialize response and components
    response = np.zeros_like(charge_axis)
    components = []
    
    # Sum over possible photoelectron counts
    for n in range(0, n_max + 1):
        # Poisson probability for n photoelectrons
        poisson_prob = poisson.pmf(n, mu_pe)
        if poisson_prob > 1e-10:  # Only consider significant contributions
            # Gaussian for n photoelectrons
            mean_n = n * Q_pe
            sigma_n = np.sqrt(n) * sigma_SPE if n > 0 else sigma_SPE
            
            # Avoid division by zero for n=0
            if sigma_n < 1e-10:
                sigma_n = sigma_SPE
                
            gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                      np.exp(-0.5 * ((charge_axis - mean_n) / sigma_n) ** 2)
            
            component = poisson_prob * gaussian
            response += component
            components.append((n, mean_n, sigma_n, component))
    
    return charge_axis, response, components

# Parameters
light_yield = 10  # photoelectrons/MeV
Q_pe = 1.0        # charge per photoelectron
sigma_SPE = 0.15  # standard deviation of SPE peak

mean_SPE = 110.342
sigma2_SPE = 10.0004

mean_pedastal = 85.5295
sigma_pedastal = 0.48595


# Create figure with subplots
fig = plt.figure(figsize=(15, 10))
gs = gridspec.GridSpec(2, 2)

# Plot 1: Different energy depositions
ax1 = plt.subplot(gs[0, 0])
PEs = [1, 5, 10, 20, 50, 100]
colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown']

mean_SPE = 110.342
sigma_SPE= 10.0004
mean_ped= 85.5295
sigma_ped= 0.48595

pedastal_idx = int(mean_ped)
PDE  = 0.2316
maxE = (4001-pedastal_idx-mean_ped)/(mean_SPE-mean_ped)/PDE

pe_axis = np.linspace(0, maxE, 4001-pedastal_idx)
print(pe_axis[1]-pe_axis[0])

for true_pe, color in zip(PEs, colors):
    charge_axis, response, _ = pmt_response(true_pe, pe_axis)
    ax1.plot(charge_axis, response, color=color, linewidth=2, 
             label=f'μ={true_pe:.1f} p.e.')

ax1.set_xlabel('PE / p.e.')
ax1.set_ylabel('Probability Density')
ax1.set_title('PMT Response for Different PE Depositions')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Individual photoelectron peaks for a specific energy
ax2 = plt.subplot(gs[0, 1])
PE = 100  # Mean number of photoelectrons
charge_axis, total_response, components = pmt_response(PE, pe_axis)

# Plot individual components
for n, mean_n, sigma_n, component in components:
    if component.max() > 0.001:  # Only plot significant components
        ax2.plot(charge_axis, component, '--', alpha=0.7, 
                label=f'N={n} (μ={mean_n:.1f}, σ={sigma_n:.3f})')

# Plot total response
ax2.plot(charge_axis, total_response, 'k-', linewidth=3, label='Total Response')

ax2.set_xlabel('Charge (arb. units)')
ax2.set_ylabel('Probability Density')
ax2.set_title(f'Individual Photoelectron Peaks ({PE} pe)')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

# Plot 3: Evolution of peak widths
ax3 = plt.subplot(gs[1, 0])
n_values = np.arange(1, 11)
theoretical_widths = np.sqrt(n_values) * sigma_SPE

ax3.plot(n_values, theoretical_widths, 'bo-', linewidth=2, markersize=6, 
         label='Theoretical: σ = √N × σ_SPE')
ax3.plot(n_values, np.sqrt(n_values) * 1.0, 'r--', alpha=0.7, 
         label='Reference: √N (normalized)')

ax3.set_xlabel('Number of Photoelectrons (N)')
ax3.set_ylabel('Peak Width σ')
ax3.set_title('Evolution of Peak Width with Photoelectron Count')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Low light level vs high light level
ax4 = plt.subplot(gs[1, 1])
# Low light
charge_axis_low, response_low, _ = pmt_response(0.3, pe_axis)
# High light  
charge_axis_high, response_high, _ = pmt_response(8.0, pe_axis)

ax4.plot(charge_axis_low, response_low, 'b-', linewidth=2, 
         label=f'Low light (E=0.3 MeV, μ={0.3*light_yield} p.e.)')
ax4.plot(charge_axis_high, response_high, 'r-', linewidth=2, 
         label=f'High light (E=8.0 MeV, μ={8.0*light_yield} p.e.)')

ax4.set_xlabel('Charge (arb. units)')
ax4.set_ylabel('Probability Density')
ax4.set_title('Low vs High Light Level Response')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Additional plot: Show the √N scaling more clearly
plt.figure(figsize=(10, 6))

n_range = np.arange(1, 21)
widths = np.sqrt(n_range) * sigma_SPE
linear = n_range * sigma_SPE  # What it would be if it scaled linearly

plt.plot(n_range, widths, 'bo-', linewidth=2, markersize=6, 
         label='Actual: σ = √N × σ_SPE')
plt.plot(n_range, linear, 'r--', linewidth=2, 
         label='Linear (for comparison): σ = N × σ_SPE')
plt.plot(n_range, np.sqrt(n_range), 'g:', alpha=0.7, 
         label='√N scaling (normalized)')

plt.xlabel('Number of Photoelectrons (N)')
plt.ylabel('Peak Width σ')
plt.title('Demonstration of √N Scaling in Peak Widths')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

# Print some key information
print("SiPM Response Function Parameters:")
print(f"Light yield: {light_yield} photoelectrons/MeV")
print(f"Charge per photoelectron (Q_pe): {Q_pe}")
print(f"SPE peak width (σ_SPE): {sigma_SPE}")
print("\nPeak width scaling:")
for n in [1, 2, 3, 4, 5, 10]:
    print(f"N={n}: σ = {np.sqrt(n):.3f} × σ_SPE = {np.sqrt(n) * sigma_SPE:.3f}")