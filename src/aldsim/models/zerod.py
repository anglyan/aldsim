from ...constants import kb
from .base import IdealDoseModel

import numpy as np

class ZeroD(IdealDoseModel):

    def __init__(self, chem, **kwargs):
        super().__init__(chem, kwargs['T'], kwargs['p'])

    def saturation_curve(self):
        print(kb)
        nu = 0.25*self.chem.site_area*self.vth*self.p/(kb*self.T)*self.chem.beta0
        t0 = 1/nu
        dt = 0.01*t0
        t_arr = np.arange(0, 5*t0, dt)
        cov_arr = 1-np.exp(-t_arr/t0)
        return t_arr, cov_arr

