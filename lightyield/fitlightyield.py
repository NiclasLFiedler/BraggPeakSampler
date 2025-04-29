import numpy as np
import matplotlib.pyplot as plt

# Replace with your actual file path
file_path = 'data/bpsCrystal34.his.txt'

# Load data (assumes two columns: x and y)
data = np.loadtxt(file_path)

# Split into x and y
x = data[:, 0]
y = data[:, 1]

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(x, y, marker='o', linestyle='-')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Data from TXT file')
plt.grid(True)
plt.tight_layout()
plt.show()
