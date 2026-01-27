#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim.core.ideal.particle.batch import PlugFlowMixed

model = PlugFlowMixed(1)
print(model.Da)

x, cov = model.saturation_curve()

import matplotlib.pyplot as pt

pt.plot(x, cov)
pt.show()
