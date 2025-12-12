#Copyright © 2024, UChicago Argonne, LLC

from .rotatingdrum import RotatingDrum, calc_coverage

class SpatialWellMixed(RotatingDrum):
    """Model for continuous particle coating under well stirred approximations.

    Implementation of a non-dimensional model for particle coating
    by atomic layer deposition for moving particles under stratified
    mixing (homogeneous mixing only on the plane perpendicular to
    the direction of movement). Precursor transport is modeled
    using the well stirred approximation.

    The model assumes a first-order irreversible Langmuir kinetics
    with the sticking probability value contained in the Damkohler
    number.

    The normalized time in the model refers to the normalized residence
    time of particles in the reactor.

    This model is formally equivalent to a batch particle coating under
    the well stirred approximation in which the normalized residence
    time is replaced by the normalized dose time.

    Args:
        Da (float) : Damkohler number

    """

def saturation_curve(Da, tmax=5, dt= 0.01):
    m = SpatialWellMixed(Da)
    return m.saturation_curve(tmax, dt)

