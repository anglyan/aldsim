#Copyright © 2024, UChicago Argonne, LLC

import numpy as np

class SpatialPlugFlow:
    """Plug flow model for particle coating using spatial ALD

    Implementation of a non-dimensional model for particle coating
    by atomic layer deposition for moving particles under stratified
    mixing (homogeneous mixing only on the plane perpendicular to
    the direction of movement). Precursor transport is modeled
    using the plug flow approximation, with both precursor and
    particles moving along the same direction.

    The model assumes a first-order irreversible Langmuir kinetics
    with the sticking probability value contained in the Damkohler
    number.

    The normalized time in the model refers to the normalized residence
    time of particles in the reactor.

    Args:
        Da (float) : Damkohler number

    """
    def __init__(self, Da):
        self.Da = Da

    def calc_coverage(self, Da=None, t=1):
        """Calculates the surface coverage

        Calculates the surface coverage for a given normalized
        dose time and Damkohler number.

        Args:
            t (float, optional): the normalized dose time
            Da (float, optional): the Damkohler number. If provided,
                it overrides the current value.
        
        Returns:
            Surface coverage
        
        """
        if Da is None:
            Da = self.Da
        else:
            self.Da = Da
        return calc_coverage(Da, t)
   
    def run(self, tmax=5, dt=0.01):
        """Runs the simulation for a

        Runs the model for a range of residence time values

        Args:
            tmax (float, optional): largest normalized residence time.
            dt (float, optional): time step value.
        
        Returns:
            A tuple of residence time, surface coverage, precursor utilization arrays
        
        """
        t = np.arange(0, tmax, dt)
        c = np.array([calc_coverage(self.Da, ti) for ti in t])
        prec = np.array([calc_precursor(self.Da, ti) for ti in t])
        return t, c, prec
    
    def saturation_curve(self, tmax=5, dt=0.01):
        """Calculates the saturation curve of the ALD process

        Calculates the saturation curve using either a default or
        user-defined set of residence times.

        Args:
            tmax (float, optional): largest normalized residence time.
            dt (float, optional): time step value.
        
        Returns:
            A tuple of residence time, surface coverage arrays
        
        """
        t = np.arange(0, tmax, dt)
        c = np.array([calc_coverage(self.Da, ti) for ti in t])
        return t, c


def calc_coverage(Da, t):
    if t == 1:
        return Da/(1+Da)
    else:
        x = np.exp(-Da*(1-t))
        return 1-(1-t)/(1-t*x)
    
def calc_precursor(Da, t):
    if t == 1:
        return 1-Da/(1+Da)
    else:
        x = np.exp(-Da*(1-t))
        return (1-t)*x/(1-t*x)


def saturation_curve(Da, tmax=5, dt= 0.01):
    m = SpatialPlugFlow(Da)
    return m.saturation_curve(tmax, dt)

