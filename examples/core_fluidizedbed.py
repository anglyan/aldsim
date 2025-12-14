"""Example using FluidizedBed class from core.particle module

This example demonstrates how to:
1. Create a FluidizedBed instance with a Damkohler number of 1
2. Calculate the saturation curve
3. Plot the saturation curve using matplotlib
"""

import matplotlib.pyplot as plt
from aldsim.core.particle import FluidizedBed

# Create a FluidizedBed instance with Damkohler number = 1
fb = FluidizedBed(Da=1.0)

# Calculate the saturation curve
t, coverage = fb.saturation_curve(tmax=5, dt=0.01)

# Plot the saturation curve
plt.figure(figsize=(8, 6))
plt.plot(t, coverage, 'b-', linewidth=2, label='Da = 1.0')
plt.xlabel('Normalized Dose Time', fontsize=12)
plt.ylabel('Surface Coverage', fontsize=12)
plt.title('Saturation Curve for Fluidized Bed ALD', fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
