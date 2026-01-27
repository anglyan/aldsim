#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim.chem import Precursor, ALDsoft
import matplotlib.pyplot as pt

prec = Precursor('TMA')
ald = ALDsoft(prec, 1e19, 0.01, 0.001, 0.8)
x, y = ald.saturation_curve(400, 10)
pt.plot(x,y)
pt.show()
