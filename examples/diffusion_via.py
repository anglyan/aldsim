from aldsim.core.diffusion import solve_until
import matplotlib.pyplot as pt
import numpy as np

AR = 20
N = 80
dt = 0.01

coverages, times = solve_until(AR, N, p_stick0=1e-3, p_rec0=1e-2, target_time=10, save_every=2)

x = 0.125+0.25*np.arange(80)

for cov, t in zip(coverages, times):
    pt.plot(x, cov, label="t={}".format(t))
pt.legend()
pt.xlabel("Normalized depth")
pt.ylabel("Surface coverage")
pt.tight_layout()
pt.show()
