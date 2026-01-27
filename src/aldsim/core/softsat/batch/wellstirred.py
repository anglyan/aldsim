#Copyright © 2024-Present, UChicago Argonne, LLC

import numpy as np
from aldsim.solvers import ode_solver

class WellStirred:
    """Implement """

    def __init__(self, D1, D2, f1, f2):
        self.D1 = D1
        self.D2 = D2
        self.f1 = f1
        self.f2 = f2
        self.a = D2/D1

    def _f(self, t, y):
        """Provides the gradient for the fraction of available sites
        of the first reaction pathway"""

        dec = 1+self.f1*self.D1*y+self.f2*self.D2*np.float_power(y, self.a)
        return - self.D1*y/dec

    def saturation_curve(self, tmax=5, dt=0.01):
        t, c, _, _, _ = self.run(tmax=tmax, dt=dt)
        return t, c

    def run(self, tmax=5, dt=0.01):
        out = ode_solver(self._f, [1], tmax, t_eval=np.arange(0,tmax,dt))
        y1 = out.y[0,:]
        y2 = np.float_power(y1, self.a)
        cov1 = 1-y1
        cov2 = 1-y2
        cov = self.f1*cov1 + self.f2*cov2
        x = 1/(1+self.f1*y1 + self.f2*y2)
        return out.t, cov, x, cov1, cov2
    
    def calc_coverage(self, t):
        raise NotImplementedError
    
