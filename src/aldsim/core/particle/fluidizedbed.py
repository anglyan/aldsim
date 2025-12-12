#Copyright © 2024, UChicago Argonne, LLC

import numpy as np

class FluidizedBed:
    """Model for batch fluidized bed reactor under plug flow approximations.

    Implementation of a non-dimensional model for particle coating
    by atomic layer deposition. It assumes a well mixed approximation for
    particle mixing and plug flow approximation for precursor transport.

    The model assumes a first-order irreversible Langmuir kinetics
    with the sticking probability value contained in the Damkohler
    number.

    Args:
        Da (float) : Damkohler number

    """

    model_kwd = ["dose", "nondim"]

    def __init__(self, Da=None):
        self.Da = Da

    def calc_coverage(self, t=1, Da=None):
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
        if Da is not None:
            self.Da = Da
        else:
            Da = self.Da
        return calc_coverage(Da, t)
    
    def saturation_curve(self, tmax=5, dt= 0.01):
        """Calculates the saturation curve of the ALD process

        Calculates the saturation using either a default or
        user-defined set of times.

        Args:
            tmax (float, optional): largest normalized dose time.
            dt (float, optional): time step value.
        
        Returns:
            A tuple of time, surface coverage arrays
        
        """
        t = np.arange(0, tmax, dt)
        Da = self.Da
        c = calc_coverage(Da, t)
        return t, c
    
    def run(self, tmax=5, dt=0.01):
        """Runs the simulation for a given or predefined amount of time

        Runs the model for a predefined or user-provided time

        Args:
            tmax (float, optional): largest normalized dose time.
            dt (float, optional): time step value.
        
        Returns:
            A tuple of time, surface coverage, precursor utilization arrays
        
        """
        t = np.arange(0, tmax, dt)
        Da = self.Da
        y = 1/Da*np.log(1+(np.exp(Da)-1)*np.exp(-Da*t))
        c = 1-y
        x = np.exp(-Da*y)
        return t, c, x


def calc_coverage(Da, t):
    """Analytical expression of the surface coverage for the PlugFlowMixed model"""
    return 1 - 1/Da*np.log(1+(np.exp(Da)-1)*np.exp(-Da*t))