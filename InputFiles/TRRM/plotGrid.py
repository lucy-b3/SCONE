import json
import numpy as np
import matplotlib.pyplot as plt

# Load the JSON file
with open('slab2DTH.json', 'r') as f:
    data = json.load(f)

# Extract the flux data
flux_data = np.array(data["flux1G"]["flux1G"])
flux_data = flux_data.reshape(50, 20, 2)  # Reshape to (50, 20, 2)

# Separate flux and uncertainty
flux = flux_data[:, :, 0]
uncertainty = flux_data[:, :, 1]

# Transpose to fix orientation
flux = flux.T
uncertainty = uncertainty.T

# Extract the boundaries
x_bounds = np.array(data["flux1G"]["XBounds"])
y_bounds = np.array(data["flux1G"]["YBounds"])

# Calculate the extent from the boundaries
x_min, x_max = x_bounds[0, 0], x_bounds[-1, 1]
y_min, y_max = y_bounds[0, 0], y_bounds[-1, 1]

# Plot the flux with correct boundaries
plt.imshow(flux, cmap='viridis', extent=[x_min, x_max, y_min, y_max], origin='lower', interpolation='nearest', aspect='equal')
plt.colorbar(label='Flux')
plt.title('Flux Distribution')
plt.xlabel('X [cm]')
plt.ylabel('Y [cm]')
plt.savefig('flux_distribution.png')  # Save the plot
plt.show()

# Plot the uncertainty with correct boundaries
plt.imshow(uncertainty, cmap='plasma', extent=[x_min, x_max, y_min, y_max], origin='lower', interpolation='nearest', aspect='equal')
plt.colorbar(label='Uncertainty')
plt.title('Uncertainty Distribution')
plt.xlabel('X [cm]')
plt.ylabel('Y [cm]')
plt.savefig('uncertainty_distribution.png')  # Save the plot
plt.show()

