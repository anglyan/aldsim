#Copyright © 2024-Present, UChicago Argonne, LLC

import numpy as np
from scipy.integrate import solve_ivp

class PlugFlowMixedNonD:
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

        dec = self.f1*self.D1*y+self.f2*self.D2*np.float_power(y, self.a)
        return - self.D1*y*(1-np.exp(-dec))/dec

    def saturation_curve(self, tmax=5):

        out = solve_ivp(self._f, [0,tmax], [1], t_eval=np.arange(0,tmax,0.01))
        y1 = out.y[0,:]
        y2 = np.float_power(y1, self.a)
        cov = self.f1*(1-y1) + self.f2*(1-y2)
        return out.t, cov

    def run(self, tmax, dt):
        raise NotImplementedError
    
    def calc_coverage(self, t):
        raise NotImplementedError
    
