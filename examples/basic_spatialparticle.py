#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim import Precursor, ALDchem
from aldsim.models import SpatialWellMixed

import matplotlib.pyplot as pt

prec = Precursor(mass=150.0)
nsites = 1e19
p_stick1 = 1e-3
chem = ALDchem(prec, nsites, p_stick1)

# Create the SpatialWellMixed model for flat surface coating
# p: precursor partial pressure (Pa), p0: carrier gas pressure (Pa)
# T: temperature (K), flow: gas flow (sccm)
# L: zone length (m), W: zone width (m)
model = SpatialWellMixed(chem, p=13.2, p0=100, T=500,
                         flow=60, L=0.02, W=0.1)

# Generate saturation curve as function of web velocity
u, theta, x = model.run(umax=10)

# Plot the saturation curve
pt.semilogx(u, theta)
pt.xlabel("Web velocity, m/s")
pt.ylabel("Fractional surface coverage")
pt.savefig("example_spatialwellmixed.png", dpi=300)
pt.show()

