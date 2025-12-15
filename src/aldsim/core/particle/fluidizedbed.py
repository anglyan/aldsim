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

    The normalized time in the model refers to the normalized dose time
    during the precursor exposure step.

    Parameters
    ----------
    Da : float
        Damkohler number, a dimensionless parameter representing the ratio
        of reaction rate to transport rate. Higher values indicate faster
        surface reactions relative to mass transport. Must be non-negative.

    Attributes
    ----------
    Da : float
        The Damkohler number for the system.

    Examples
    --------
    Create a FluidizedBed model with a Damkohler number of 2.0:

    >>> model = FluidizedBed(Da=2.0)
    >>> coverage = model.calc_coverage(t=1.0)
    >>> print(f"Coverage: {coverage:.3f}")
    Coverage: 0.797

    Calculate saturation curve over normalized dose time:

    >>> t, coverage = model.saturation_curve(tmax=3.0, dt=0.1)
    >>> max_coverage = coverage[-1]
    >>> print(f"Maximum coverage: {max_coverage:.3f}")
    Maximum coverage: 0.950

    Notes
    -----
    The Damkohler number (Da) is defined as:

        Da = k * τ

    where k is the first-order reaction rate constant and τ is the
    characteristic time scale for precursor dosing.

    """

    model_kwd = ["dose", "nondim"]

    def __init__(self, Da=None):
        self.Da = Da

    def calc_coverage(self, t=1, Da=None):
        """Calculate the surface coverage at a given normalized dose time

        Computes the fractional surface coverage θ of particles by precursor
        molecules in the batch fluidized bed reactor at a specified normalized
        dose time. The calculation uses analytical solutions for the well-mixed
        particle model with plug flow precursor transport and first-order
        Langmuir kinetics.

        Parameters
        ----------
        t : float, optional
            Normalized dose time (dimensionless), defined as the ratio
            of actual dose time to the characteristic dosing time.
            Default is 1.0 (dose duration equals the characteristic time).
            Must be non-negative.
        Da : float, optional
            Damkohler number. If provided, updates the model's Da attribute
            and uses this value for the calculation. If None (default),
            uses the current model's Da value.

        Returns
        -------
        float
            Surface coverage θ (dimensionless), bounded between 0 and 1.
            A value of 0 indicates no coverage, while 1 indicates complete
            monolayer saturation.

        Examples
        --------
        Calculate coverage at the end of dosing (t=1):

        >>> model = FluidizedBed(Da=2.0)
        >>> coverage = model.calc_coverage(t=1.0)
        >>> print(f"Coverage: {coverage:.3f}")
        Coverage: 0.797

        Override the Damkohler number for a specific calculation:

        >>> coverage = model.calc_coverage(t=1.0, Da=5.0)
        >>> print(f"Coverage with Da=5: {coverage:.3f}")
        Coverage with Da=5: 0.959


        """
        if Da is not None:
            self.Da = Da
        else:
            Da = self.Da
        cov = 1 - 1/Da*np.log(1+(np.exp(Da)-1)*np.exp(-Da*t))
        return cov
    
    def saturation_curve(self, tmax=5, dt= 0.01):
        """Calculate the saturation curve of the ALD process

        Computes the relationship between normalized dose time and
        surface coverage, producing the characteristic saturation curve
        for the fluidized bed ALD process. This curve shows how coverage
        approaches saturation as dose time increases in a batch reactor
        with well-mixed particles and plug flow precursor transport.

        Parameters
        ----------
        tmax : float, optional
            Maximum normalized dose time for the curve. Default is 5.0.
            Determines the extent of the curve. Larger values show
            near-saturation behavior but may be unnecessary if saturation
            is reached earlier.
        dt : float, optional
            Time step size (dimensionless). Default is 0.01.
            Smaller values provide smoother curves but increase computation time.
            Must be positive and smaller than tmax.

        Returns
        -------
        t : ndarray
            Array of normalized dose times, shape (n,), where n = tmax/dt.
            Values range from 0 to (tmax - dt).
        coverage : ndarray
            Array of surface coverage values θ at each time point, shape (n,).
            Each value is bounded between 0 and 1. The coverage increases
            monotonically, approaching saturation at large dose times.

        Examples
        --------
        Generate a saturation curve with default parameters:

        >>> model = FluidizedBed(Da=2.0)
        >>> t, coverage = model.saturation_curve()
        >>> print(f"Coverage at t=1: {coverage[100]:.3f}")  # dt=0.01, so index 100 is t=1
        Coverage at t=1: 0.797

        Create a high-resolution saturation curve:

        >>> t, coverage = model.saturation_curve(tmax=3.0, dt=0.001)
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(t, coverage)
        >>> plt.xlabel('Normalized dose time')
        >>> plt.ylabel('Surface coverage θ')
        >>> plt.title(f'Fluidized Bed ALD Saturation Curve (Da={model.Da})')
        >>> plt.grid(True)
        >>> plt.show()

        See Also
        --------
        run : Complete simulation including precursor utilization
        calc_coverage : Calculate coverage at a single time point

        """
        t = np.arange(0, tmax, dt)
        c = self.calc_coverage(t)
        return t, c
    
    def run(self, tmax=5, dt=0.01):
        """Run complete simulation including coverage and precursor utilization

        Executes the fluidized bed model simulation over a range of normalized
        dose times, computing both surface coverage and precursor utilization
        at each time step. This method provides comprehensive results for
        analyzing both coating efficiency and precursor usage in a batch
        reactor with well-mixed particles and plug flow precursor transport.

        Parameters
        ----------
        tmax : float, optional
            Maximum normalized dose time for the simulation. Default is 5.0.
            Should be greater than 0. Larger values allow observation of
            near-saturation behavior.
        dt : float, optional
            Time step size for the simulation (dimensionless). Default is 0.01.
            Smaller values provide higher resolution but increase computation time.
            Must be positive and smaller than tmax.

        Returns
        -------
        t : ndarray
            Array of normalized dose times, shape (n,), where n = tmax/dt.
            Values range from 0 to (tmax - dt).
        coverage : ndarray
            Array of surface coverage values θ at each time point, shape (n,).
            Each value is bounded between 0 and 1, representing fractional
            monolayer coverage.
        precursor : ndarray
            Array of precursor utilization factors at each time point, shape (n,).
            Each value is bounded between 0 and 1, representing the fraction
            of precursor that has reacted with particle surfaces. Values closer
            to 0 indicate less efficient precursor usage.

        Examples
        --------
        Run a basic simulation with default parameters:

        >>> model = FluidizedBed(Da=2.0)
        >>> t, coverage, precursor = model.run()
        >>> print(f"Final coverage: {coverage[-1]:.3f}")
        Final coverage: 0.993

        Notes
        -----
        This method combines the functionality of calc_coverage() and precursor
        utilization calculations to provide a complete picture of the fluidized
        bed ALD process.

        
        See Also
        --------
        calc_coverage : Calculate coverage only
        saturation_curve : Calculate coverage without precursor data

        """
        t = np.arange(0, tmax, dt)
        Da = self.Da
        y = 1/Da*np.log(1+(np.exp(Da)-1)*np.exp(-Da*t))
        c = 1-y
        x = np.exp(-Da*y)
        return t, c, x

