import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.signal import convolve
import matplotlib.gridspec as gridspec

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
energies = [0.5, 1.0, 2.0, 4.0]  # MeV
colors = ['blue', 'green', 'red', 'purple']

for energy, color in zip(energies, colors):
    charge_axis, response, _ = sipm_response(energy, light_yield, Q_pe, sigma_SPE)
    ax1.plot(charge_axis, response, color=color, linewidth=2, 
             label=f'{energy} MeV (μ={energy*light_yield:.1f} p.e.)')

ax1.set_xlabel('Charge (arb. units)')
ax1.set_ylabel('Probability Density')
ax1.set_title('SiPM Response for Different Energy Depositions')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Individual photoelectron peaks for a specific energy
ax2 = plt.subplot(gs[0, 1])
energy = 2.0  # MeV
charge_axis, total_response, components = sipm_response(energy, light_yield, Q_pe, sigma_SPE)

# Plot individual components
for n, mean_n, sigma_n, component in components:
    if component.max() > 0.001:  # Only plot significant components
        ax2.plot(charge_axis, component, '--', alpha=0.7, 
                label=f'N={n} (μ={mean_n:.1f}, σ={sigma_n:.3f})')

# Plot total response
ax2.plot(charge_axis, total_response, 'k-', linewidth=3, label='Total Response')

ax2.set_xlabel('Charge (arb. units)')
ax2.set_ylabel('Probability Density')
ax2.set_title(f'Individual Photoelectron Peaks (E={energy} MeV)')
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
charge_axis_low, response_low, _ = sipm_response(0.3, light_yield, Q_pe, sigma_SPE)
# High light  
charge_axis_high, response_high, _ = sipm_response(8.0, light_yield, Q_pe, sigma_SPE)

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