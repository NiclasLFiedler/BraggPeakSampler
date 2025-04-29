import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import linregress
matplotlib.use('TkAgg')  # or 'Qt5Agg'

# Sample data points
x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # X values
y = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11])  # Y values

ch_Ba133 = [762.5, 1507, 2093, 2594, 2842, 3345, 3605, 4110]
ch_Eu152 = [1144, 2294, 3232, 3449, 3862, 4173, 7331, 8167, 9074, 10226, 10471, 11417, 12232, 13255]
ch_Ra226 = [1745, 2270, 2771, 3308, 5738, 6270, 7241, 8802, 10558, 11665, 12979, 14216]

e_Ba133 = [81, 160.6, 223.2, 276.4, 302.8, 356, 383.8, 437]
e_Eu152 = [121.8, 244.7, 344.3, 367.8, 411.1, 444, 778.9, 867.4, 964.1, 1085.9, 1112.1, 1212.9, 1299, 1408]
e_Ra226 = [186, 241.9, 295.2, 351.9, 609.3, 665.5, 768.4, 934.1, 1120.3, 1238.1, 1377.7, 1509]

# Combine all data sets into one
ch_all = np.concatenate([ch_Ba133, ch_Eu152, ch_Ra226])
e_all = np.concatenate([e_Ba133, e_Eu152, e_Ra226])

# Perform linear regression on the combined dataset
slope, intercept, r_value, p_value, std_err = linregress(ch_all, e_all)
r_squared = r_value**2
fit_line = slope * ch_all + intercept

# Create the plot
plt.figure(figsize=(8, 7))

# Plot all data points
plt.scatter(ch_Ba133, e_Ba133, color='blue', label='Ba-133 Data Points', zorder=5)
plt.scatter(ch_Eu152, e_Eu152, color='orange', label='Eu-152 Data Points', zorder=5)
plt.scatter(ch_Ra226, e_Ra226, color='green', label='Ra-226 Data Points', zorder=5)

# Plot the linear fit
plt.plot(ch_all, fit_line, color='red', label=f'Linear fit: $R^2={r_squared:.8f}$ \nSlope: {slope:.4e}\nIntercept:{intercept:.4e}', linewidth=2, zorder=10)

# Enhance plot
plt.title('Linearity fit of CAEN V1730S digitizer', fontsize=16)
plt.xlabel('ADC Channel', fontsize=14)
plt.ylabel('Energy / keV', fontsize=14)
plt.grid(True)
plt.legend()
plt.tight_layout()

# Show the plot
plt.savefig("linearity.pdf", dpi=300, bbox_inches='tight', format="pdf")
plt.show()

# Print R² value
print(f"R² value: {r_squared:.9f}")
for ch in ch_Ra226:
    E = slope * ch + intercept
    print(f"{E:.3f}")