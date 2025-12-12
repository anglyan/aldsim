from aldsim import Precursor, ALDideal, aldmodel

import matplotlib.pyplot as pt

prec = Precursor(mass=150.0)
nsites = 1e19
sprob0 = 1e-3
chem = ALDideal(prec, nsites, sprob0, dm=1.0)
model = aldmodel(chem, 'fluidizedbed', p=0.1*1e5/760, p0=1e2, T=500,
                 S=1e1, flow=60)

t, theta = model.saturation_curve()

pt.plot(t, theta)
pt.show()

