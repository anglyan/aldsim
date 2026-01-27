#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim import Precursor, ALDideal, aldmodel

import matplotlib.pyplot as pt

prec = Precursor(mass=150.0)
nsites = 1e19
beta0 = 1e-3
chem = ALDideal(prec, nsites, beta0, dm=1.0)
model = aldmodel(chem, 'zeroD', p=0.1*1e5/760, T=500)
t, theta = model.saturation_curve()

pt.plot(t, theta)
pt.show()



