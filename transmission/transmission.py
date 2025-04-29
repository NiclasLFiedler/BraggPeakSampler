import matplotlib.pyplot as plt
import pandas as pd

# Load the data
file_path = "PWO_transmission2.dat"  # Replace with your file name
data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=["Wavelength / nm", "Transmittance / %"], decimal=",")

# Plot the data
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(10, 6))
plt.plot(data["Wavelength / nm"], data["Transmittance / %"], color='blue', linewidth=1.5)

# Customize plot appearance
plt.xlabel("Wavelength / nm")
plt.ylabel("Transmittance / %")
plt.grid(True)

# Add a legend and show the plot
plt.legend(["Transmittance"], loc="upper left")
plt.tight_layout()
plt.savefig("pwo_transmission.pdf", format="pdf")
plt.show()
