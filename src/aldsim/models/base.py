#Copyright © 2024-Present, UChicago Argonne, LLC

class DoseModel:

    def __init__(self, chem, T, p):
        self.chem = chem
        self.p = p
        self.T = T

    @property
    def T(self):
        return self._T
 
    @T.setter
    def T(self, value):
        self._vth = self.chem.vth(value)
        self._T = value

    @property
    def vth(self):
        return self._vth
    
    @property
    def mass(self):
        return self.prec.mass
   
    @property
    def p(self):
        return self._p
 
    @p.setter
    def p(self, value):
        self._p = value

    def run(self, time):
        raise NotImplementedError("Models require a run method")
    
    