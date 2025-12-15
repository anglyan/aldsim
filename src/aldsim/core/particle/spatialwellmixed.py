#Copyright © 2024, UChicago Argonne, LLC

from .rotatingdrum import RotatingDrum

class SpatialWellMixed(RotatingDrum):
    """Model for continuous particle coating under well stirred approximations.

    Implementation of a non-dimensional model for particle coating
    by atomic layer deposition for moving particles under stratified
    mixing (homogeneous mixing only on the plane perpendicular to
    the direction of movement). Precursor transport is modeled
    using the well stirred approximation. This model is applicable to
    continuous flow systems where particles move through the reactor
    while the gas phase is well stirred.

    The model assumes a first-order irreversible Langmuir kinetics
    with the sticking probability value contained in the Damkohler
    number.

    The normalized time in the model refers to the normalized residence
    time of particles in the reactor.

    This model is formally equivalent to a batch particle coating under
    the well stirred approximation (RotatingDrum) in which the normalized
    residence time is replaced by the normalized dose time. The mathematical
    formulation is identical, allowing SpatialWellMixed to inherit all
    methods from RotatingDrum.

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
    Create a SpatialWellMixed model with a Damkohler number of 2.0:

    >>> model = SpatialWellMixed(Da=2.0)
    >>> coverage = model.calc_coverage(t=1.0)
    >>> print(f"Coverage: {coverage:.3f}")
    Coverage: 0.757

    Calculate saturation curve over normalized residence time:

    >>> t, coverage = model.saturation_curve(tmax=3.0, dt=0.1)
    >>> max_coverage = coverage[-1]
    >>> print(f"Maximum coverage: {max_coverage:.3f}")
    Maximum coverage: 0.950

    See Also
    --------
    RotatingDrum : Batch reactor model with identical mathematical formulation
    SpatialPlugFlow : Continuous model with plug flow precursor transport

    """


