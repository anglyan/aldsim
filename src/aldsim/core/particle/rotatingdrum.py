#Copyright © 2024, UChicago Argonne, LLC

import numpy as np
from aldsim.solvers import ode_solver, boundedNewton_solver

class RotatingDrum:
    """Model for batch particle coating under a well mixed reactor approximation.

    Implementation of a non-dimensional model for particle coating
    by atomic layer deposition under a well mixed approximation for
    particle mixing and well stirred approximation for precursor transport.

    The model assumes a first-order irreversible Langmuir kinetics
    with the sticking probability value contained in the Damkohler
    number.

    Args:
        Da (float) : Damkohler number

    """
    def __init__(self, Da=None):
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
    
    def _f_t(self, theta, Da, tau):
        return theta - np.log(1-theta)/Da - tau

    def _fp_t(self, theta, Da, tau):
        return 1 + 1/(Da*(1-theta))

    def _f(self, t, y):
        return -y/(1/self.Da+y)

    def run(self, tmax=5, dt=0.01):
        """Runs the simulation for a given or predefined amount of time

        Runs the model for a predefined or user-provided time

        Args:
            tmax (float, optional): largest normalized dose time.
            dt (float, optional): time step value.
        
        Returns:
            A tuple of time, surface coverage, precursor utilization arrays
        
        """
        out = ode_solver(self._f, [1], tmax, t_eval=np.arange(0,tmax,dt))
        cov = 1-out.y[0,:]
        x = 1/(1+self.Da*out.y[0,:])
        return out.t, cov, x

    def saturation_curve(self, tmax=5, dt=0.01):
        """Calculates the saturation curve of the ALD process

        Calculates the saturation using either a default or
        user-defined set of times.

        Args:
            tmax (float, optional): largest normalized dose time.
            dt (float, optional): time step value.
        
        Returns:
            A tuple of time, surface coverage arrays
        
        """
        t, cov, _ = self.run(tmax, dt)
        return t, cov
        
    def saturation_curve_implicit(self, theta_max=0.9999):
        Da = self.Da
        theta = np.arange(0,theta_max,0.0001)
        tau = theta - np.log(1-theta)/Da
        return tau, theta

    def fraction_out(self, theta_max=0.999):
        Da = self.Da
        theta = np.arange(0,theta_max,0.0001)
        tau = theta - np.log(1-theta)/Da
        return tau, 1/(1+Da*(1-theta))


def calc_coverage(Da, tau):

    f = lambda t: t - np.log(1-t)/Da - tau
    fdot = lambda t: 1 + 1/(Da*(1-t))
    t = boundedNewton_solver(f, fdot)
    return t


def saturation_curve_double(Da1, Da2, f1, f2, theta_max=0.99999):
    alpha = Da2/Da1
    theta2 = np.arange(0,theta_max,0.00001)
    x2 = 1-theta2
    x1 = np.power(x2, 1/alpha)
    theta1 = 1-x1
    tau = -np.log(x1)/Da1 + f1*theta1 + f2*(1-np.power(x1, alpha))
    return tau, f1*theta1+f2*theta2


def saturation_curve(Da, tmax=5, dt= 0.01):
    m = RotatingDrum(Da)
    return m.saturation_curve(tmax, dt)
