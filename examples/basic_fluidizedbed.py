#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim import Precursor, ALDchem
from aldsim.models import FluidizedBed
import matplotlib.pyplot as pt


# Define precursor and surface chemistry
tma = Precursor('TMA', mass=144.17)
chem = ALDchem(tma, nsites=1e19, p_stick1=1e-3)

# Create the FluidizedBed model
# p: precursor partial pressure (Pa), p0: carrier gas pressure (Pa)
# T: temperature (K), S: total particle surface area (m²), flow: gas flow (sccm)
model = FluidizedBed(chem, p=13.2, p0=100, T=500, S=10, flow=60)

# Generate the saturation curve
time, coverage = model.saturation_curve()

# Plot the saturation curve
pt.plot(time, coverage)
pt.xlabel("Dose time, s")
pt.ylabel("Fractional surface coverage")
pt.savefig("example_fluidizedbed.png", dpi=300)
pt.show()

