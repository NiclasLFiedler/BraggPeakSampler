import numpy as np
import matplotlib.pyplot as plt

filename = '03-Fricke-20Gy_7514292SP_7514292SP_0030.Raw8'
width = 640                 
height = 480                

with open(filename, 'rb') as f:
    data = np.fromfile(f, dtype=np.uint8)

if data.size == width * height:
    data = data.reshape((height, width))


# Read binary data
raw_data = np.fromfile(filename, dtype=np.uint8)

# Reshape
image = raw_data.reshape((height, width))

# Display image
plt.imshow(image, cmap='gray', vmin=0, vmax=255)
plt.title("RAW8 Image")
plt.axis('off')
plt.show()
