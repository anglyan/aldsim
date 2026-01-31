#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim import Precursor, ALDchem
from aldsim.models import ZeroD
import matplotlib.pyplot as pt

# Define a precursor (TMA with molecular mass in g/mol)
tma = Precursor('TMA', mass=144.17)

# Define surface kinetics with ideal Langmuir behavior
kinetics = ALDchem(prec=tma, nsites=1e19, p_stick1=0.01)

# Create the ZeroD model at 200°C and 100 Pa
model = ZeroD(chem=kinetics, T=473.15, p=100)

# Generate the saturation curve
time, coverage = model.saturation_curve()
pt.plot(time, coverage)
pt.xlabel("Dose time, s")
pt.ylabel("Fractional surface coverage")
pt.savefig("example_zerod.png", dpi=300)
pt.show()
