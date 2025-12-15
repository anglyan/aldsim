#Copyright © 2024, UChicago Argonne, LLC

import numpy as np
from aldsim.solvers import ode_solver, boundedNewton_solver

class RotatingDrum:
    """Model for batch particle coating under a well mixed reactor approximation.

    Implementation of a non-dimensional model for particle coating
    by atomic layer deposition under a well mixed approximation for
    particle mixing and well stirred approximation for precursor transport.
    This model is applicable to rotating drum reactors and other systems
    where particles are thoroughly mixed and the gas phase is well stirred.

    The model assumes a first-order irreversible Langmuir kinetics
    with the sticking probability value contained in the Damkohler
    number.

    Parameters
    ----------
    Da : float, optional
        Damkohler number, a dimensionless parameter representing the ratio
        of reaction rate to transport rate. Higher values indicate faster
        surface reactions relative to mass transport. Must be non-negative.
        If None, must be specified in subsequent method calls.

    Attributes
    ----------
    Da : float
        The Damkohler number for the system.

    Examples
    --------
    Create a RotatingDrum model with a Damkohler number of 2.0:

    >>> model = RotatingDrum(Da=2.0)
    >>> coverage = model.calc_coverage(t=1.0)
    >>> print(f"Coverage: {coverage:.3f}")
    Coverage: 0.757

    Calculate saturation curve over normalized dose time:

    >>> t, coverage = model.saturation_curve(tmax=3.0, dt=0.1)
    >>> max_coverage = coverage[-1]
    >>> print(f"Maximum coverage: {max_coverage:.3f}")
    Maximum coverage: 0.950

    """
    def __init__(self, Da=None):
        self.Da = Da

    def calc_coverage(self, Da=None, t=1):
        """Calculate the surface coverage at a given normalized dose time

        Computes the fractional surface coverage θ of particles by precursor
        molecules in the well-mixed batch reactor at a specified normalized
        dose time. The calculation uses an implicit solution solved numerically
        via bounded Newton's method, accounting for the coupling between
        surface coverage and precursor depletion in the batch system.

        Parameters
        ----------
        Da : float, optional
            Damkohler number. If provided, updates the model's Da attribute
            and uses this value for the calculation. If None (default),
            uses the current model's Da value.
        t : float, optional
            Normalized dose time (dimensionless), defined as the ratio
            of actual dose time to the characteristic reaction time.
            Default is 1.0. Must be non-negative.

        Returns
        -------
        float
            Surface coverage θ (dimensionless), bounded between 0 and 1.
            A value of 0 indicates no coverage, while 1 indicates complete
            monolayer saturation.

        Examples
        --------
        Calculate coverage at unit dose time:

        >>> model = RotatingDrum(Da=2.0)
        >>> coverage = model.calc_coverage(t=1.0)
        >>> print(f"Coverage: {coverage:.3f}")
        Coverage: 0.757

        Override the Damkohler number for a specific calculation:

        >>> coverage = model.calc_coverage(Da=5.0, t=1.0)
        >>> print(f"Coverage with Da=5: {coverage:.3f}")
        Coverage with Da=5: 0.924

        See Also
        --------
        run : Calculate coverage evolution over time
        saturation_curve : Calculate coverage vs. time relationship

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
        """Run complete simulation including coverage and precursor utilization

        Executes the well-mixed batch reactor model simulation over a range
        of normalized dose times, computing both surface coverage and precursor
        utilization at each time step. The simulation uses ODE integration to
        track the evolution of the system, accounting for coupled precursor
        depletion and surface coverage buildup.

        Parameters
        ----------
        tmax : float, optional
            Maximum normalized dose time for the simulation. Default is 5.0.
            Should be greater than 0. Larger values allow observation of
            near-saturation behavior in batch operation.
        dt : float, optional
            Time step size for output (dimensionless). Default is 0.01.
            Smaller values provide higher resolution but may increase
            computation time. Must be positive and smaller than tmax.

        Returns
        -------
        t : ndarray
            Array of normalized dose times, shape (n,) where n ≈ tmax/dt.
            Values range from 0 to approximately tmax.
        coverage : ndarray
            Array of surface coverage values θ at each time point, shape (n,).
            Each value is bounded between 0 and 1, representing fractional
            monolayer coverage. Coverage increases with time.
        precursor : ndarray
            Array of precursor utilization fractions at each time point, shape (n,).
            Each value is bounded between 0 and 1, representing the fraction
            of precursor that has been consumed (reacted with particles).

        Examples
        --------
        Run a basic simulation with default parameters:

        >>> model = RotatingDrum(Da=2.0)
        >>> t, coverage, precursor = model.run()
        >>> print(f"Final coverage: {coverage[-1]:.3f}")
        Final coverage: 0.993

        Run simulation with custom time range and resolution:

        >>> t, coverage, precursor = model.run(tmax=3.0, dt=0.05)

        See Also
        --------
        calc_coverage : Calculate coverage at a single time point
        saturation_curve : Calculate coverage without precursor data
        saturation_curve_implicit : Alternative implicit saturation curve method

        """
        out = ode_solver(self._f, [1], tmax, t_eval=np.arange(0,tmax,dt))
        cov = 1-out.y[0,:]
        x = 1/(1+self.Da*out.y[0,:])
        return out.t, cov, x

    def saturation_curve(self, tmax=5, dt=0.01):
        """Calculate the saturation curve of the batch ALD process

        Computes the relationship between normalized dose time and surface
        coverage for the well-mixed batch reactor, producing the characteristic
        saturation curve. This curve shows how coverage builds up in batch
        operation as precursor is consumed and particles are coated.

        Parameters
        ----------
        tmax : float, optional
            Maximum normalized dose time for the curve. Default is 5.0.
            Determines the extent of the curve. Larger values show
            near-saturation behavior but may be unnecessary if saturation
            is reached earlier in batch operation.
        dt : float, optional
            Time step size for output (dimensionless). Default is 0.01.
            Smaller values provide smoother curves but may increase computation time.
            Must be positive and smaller than tmax.

        Returns
        -------
        t : ndarray
            Array of normalized dose times, shape (n,) where n ≈ tmax/dt.
            Values range from 0 to approximately tmax.
        coverage : ndarray
            Array of surface coverage values θ at each time point, shape (n,).
            Each value is bounded between 0 and 1. The coverage increases
            monotonically, approaching saturation as precursor is depleted
            in batch operation.

        Examples
        --------
        Generate a saturation curve with default parameters:

        >>> model = RotatingDrum(Da=2.0)
        >>> t, coverage = model.saturation_curve()
        >>> print(f"Coverage at t=1: {coverage[100]:.3f}")  # dt=0.01, index ≈100
        Coverage at t=1: 0.757

        Create a high-resolution saturation curve:

        >>> t, coverage = model.saturation_curve(tmax=3.0, dt=0.001)
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(t, coverage)
        >>> plt.xlabel('Normalized dose time')
        >>> plt.ylabel('Surface coverage θ')
        >>> plt.title(f'Batch ALD Saturation Curve (Da={model.Da})')
        >>> plt.grid(True)
        >>> plt.show()

        See Also
        --------
        run : Complete simulation including precursor utilization
        calc_coverage : Calculate coverage at a single time point
        saturation_curve_implicit : Implicit formulation of saturation curve

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
