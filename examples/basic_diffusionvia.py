#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim import Precursor, ALDchem
from aldsim.models import KnudsenVia
import matplotlib.pyplot as pt

# Define a precursor (TMA with molecular mass in g/mol)
tma = Precursor('TMA', mass=144.17)

# Define surface kinetics with ideal Langmuir behavior
kinetics = ALDchem(prec=tma, nsites=1e19, p_stick1=0.01)

# Create the ZeroD model at 200°C and 10 Pa
model = KnudsenVia(chem=kinetics, T=473.15, p=10, AR=200)

# Generate the saturation curve
time, z, coverage = model.saturation_curve()


for t, c in zip(time, coverage):
    pt.plot(z, c, label=("%4.2f s"% t))

pt.xlabel("Depth z/d")
pt.ylabel("Coverage")
pt.legend()
pt.savefig("example_circularvia.png", dpi=300)
pt.show()
