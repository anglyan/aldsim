from aldsim.core.diffusion import DiffusionVia
import matplotlib.pyplot as pt
import numpy as np

AR = 20
N = 80
dt = 0.01

dv = DiffusionVia(AR, p_stick0=1e-3)

coverage, times = dv.run_until_cov(max_cov=0.5, N=N, save_every=0.2, dt=dt)

x = 0.125+0.25*np.arange(80)

for cov, t in zip(coverage, times):
    pt.plot(x, cov, label="t={}".format(t))
pt.legend()
pt.xlabel("Normalized depth")
pt.ylabel("Surface coverage")
pt.tight_layout()
pt.show()
