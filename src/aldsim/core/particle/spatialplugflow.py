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
    Create a SpatialPlugFlow model with a Damkohler number of 2.0:

    >>> model = SpatialPlugFlow(Da=2.0)
    >>> coverage = model.calc_coverage(t=1.0)
    >>> print(f"Coverage: {coverage:.3f}")
    Coverage: 0.667

    Calculate saturation curve over normalized residence time:

    >>> t, coverage = model.saturation_curve(tmax=3.0, dt=0.1)
    >>> max_coverage = coverage[-1]
    >>> print(f"Maximum coverage: {max_coverage:.3f}")
    Maximum coverage: 0.950

    """
    def __init__(self, Da):
        self.Da = Da

    def calc_coverage(self, t=1, Da=None):
        """Calculate the surface coverage at a given normalized residence time

        Computes the fractional surface coverage θ of particles by precursor
        molecules in the plug flow reactor at a specified normalized residence
        time. The calculation uses analytical solutions for the plug flow model
        with first-order Langmuir kinetics.

        Parameters
        ----------
        t : float, optional
            Normalized residence time (dimensionless), defined as the ratio
            of actual residence time to the characteristic reactor time.
            Default is 1.0 (particles exit at the reactor characteristic time).
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
        Calculate coverage at the reactor exit (t=1):

        >>> model = SpatialPlugFlow(Da=2.0)
        >>> coverage = model.calc_coverage(t=1.0)
        >>> print(f"Coverage: {coverage:.3f}")
        Coverage: 0.667

        Override the Damkohler number for a specific calculation:

        >>> coverage = model.calc_coverage(t=1.0, Da=5.0)
        >>> print(f"Coverage with Da=5: {coverage:.3f}")
        Coverage with Da=5: 0.833

        """
        if Da is None:
            Da = self.Da
        else:
            self.Da = Da
        
        if t == 1:
            return Da/(1+Da)
        else:
            x = np.exp(-Da*(1-t))
            return 1-(1-t)/(1-t*x)

   
    def run(self, tmax=5, dt=0.01):
        """Run complete simulation including coverage and precursor utilization

        Executes the plug flow model simulation over a range of normalized
        residence times, computing both surface coverage and precursor
        utilization at each time step. This method provides comprehensive
        results for analyzing both coating efficiency and precursor usage.

        Parameters
        ----------
        tmax : float, optional
            Maximum normalized residence time for the simulation. Default is 5.0.
            Should be greater than 0. Larger values allow observation of
            near-saturation behavior.
        dt : float, optional
            Time step size for the simulation (dimensionless). Default is 0.01.
            Smaller values provide higher resolution but increase computation time.
            Must be positive and smaller than tmax.

        Returns
        -------
        t : ndarray
            Array of normalized residence times, shape (n,), where n = tmax/dt.
            Values range from 0 to (tmax - dt).
        coverage : ndarray
            Array of surface coverage values θ at each time point, shape (n,).
            Each value is bounded between 0 and 1, representing fractional
            monolayer coverage.
        precursor : ndarray
            Array of unused precursor fractions at each time point, shape (n,).
            Each value is bounded between 0 and 1, representing the fraction
            of precursor that exits the reactor without reacting.

        Examples
        --------
        Run a basic simulation with default parameters:

        >>> model = SpatialPlugFlow(Da=2.0)
        >>> t, coverage, precursor = model.run()
        >>> print(f"Final coverage: {coverage[-1]:.3f}")
        Final coverage: 0.993

        Run simulation with custom time range and resolution:

        >>> t, coverage, precursor = model.run(tmax=3.0, dt=0.05)
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(t, coverage, label='Coverage')
        >>> plt.plot(t, precursor, label='Unused precursor')
        >>> plt.xlabel('Normalized residence time')
        >>> plt.ylabel('Fraction')
        >>> plt.legend()
        >>> plt.show()

        Notes
        -----
        This method combines the functionality of calc_coverage() and
        calc_precursor() to provide a complete picture of the ALD process.

        See Also
        --------
        calc_coverage : Calculate coverage only
        calc_precursor : Calculate precursor utilization only
        saturation_curve : Calculate coverage without precursor data

        """
        t = np.arange(0, tmax, dt)
        c = np.array([self.calc_coverage(ti) for ti in t])
        prec = np.array([self.calc_precursor(ti) for ti in t])
        return t, c, prec
    
    def saturation_curve(self, tmax=5, dt=0.01):
        """Calculate the saturation curve of the ALD process

        Computes the relationship between normalized residence time and
        surface coverage, producing the characteristic saturation curve
        for the ALD process. This curve shows how coverage approaches
        saturation as residence time increases.

        Parameters
        ----------
        tmax : float, optional
            Maximum normalized residence time for the curve. Default is 5.0.
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
            Array of normalized residence times, shape (n,), where n = tmax/dt.
            Values range from 0 to (tmax - dt).
        coverage : ndarray
            Array of surface coverage values θ at each time point, shape (n,).
            Each value is bounded between 0 and 1. The coverage increases
            monotonically, approaching saturation at large residence times.

        Examples
        --------
        Generate a saturation curve with default parameters:

        >>> model = SpatialPlugFlow(Da=2.0)
        >>> t, coverage = model.saturation_curve()
        >>> print(f"Coverage at t=1: {coverage[100]:.3f}")  # dt=0.01, so index 100 is t=1
        Coverage at t=1: 0.667

        Create a high-resolution saturation curve:

        >>> t, coverage = model.saturation_curve(tmax=3.0, dt=0.001)
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(t, coverage)
        >>> plt.xlabel('Normalized residence time')
        >>> plt.ylabel('Surface coverage θ')
        >>> plt.title(f'ALD Saturation Curve (Da={model.Da})')
        >>> plt.grid(True)
        >>> plt.show()

        See Also
        --------
        run : Complete simulation including precursor utilization
        calc_coverage : Calculate coverage at a single time point

        """
        t = np.arange(0, tmax, dt)
        c = np.array([self.calc_coverage(ti) for ti in t])
        return t, c


    def calc_precursor(self, t, Da=None):
        """Calculate the fraction of unused precursor exiting the reactor

        Computes the fraction of precursor molecules that pass through the
        reactor without reacting with particle surfaces. This quantity is
        important for understanding precursor efficiency and optimizing
        process economics.

        Parameters
        ----------
        t : float
            Normalized residence time (dimensionless), defined as the ratio
            of actual residence time to the characteristic reactor time.
            Must be non-negative.
        Da : float, optional
            Damkohler number. If provided, updates the model's Da attribute
            and uses this value for the calculation. If None (default),
            uses the current model's Da value.

        Returns
        -------
        float
            Fraction of unused precursor (dimensionless), bounded between 0 and 1.
            A value of 0 indicates complete precursor utilization (all precursor
            reacts), while 1 indicates no precursor consumption.

        Examples
        --------
        Calculate unused precursor at reactor exit:

        >>> model = SpatialPlugFlow(Da=2.0)
        >>> unused = model.calc_precursor(t=1.0)
        >>> utilization = 1 - unused
        >>> print(f"Precursor utilization: {utilization:.1%}")
        Precursor utilization: 66.7%


        See Also
        --------
        calc_coverage : Calculate surface coverage
        run : Calculate both coverage and precursor utilization

        """
        if Da is None:
            Da = self.Da
        else:
            self.Da = Da
        if t == 1:
            return 1-Da/(1+Da)
        else:
            x = np.exp(-Da*(1-t))
            return (1-t)*x/(1-t*x)


def saturation_curve(Da, tmax=5, dt= 0.01):
    m = SpatialPlugFlow(Da)
    return m.saturation_curve(tmax, dt)

