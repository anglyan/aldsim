#Copyright © 2024-Present, UChicago Argonne, LLC

from ..constants import kb
from .base import DoseModel

import numpy as np

class ZeroD(DoseModel):
    """
    Zero-dimensional model for the time evolution of a self-limited surface process.

    This class models the transient behavior of a self-limited reaction (such as ALD)
    at constant precursor pressure and temperature. The model assumes spatially uniform
    conditions (no transport limitations), making it suitable for describing ideal
    surface kinetics or well-mixed reactor conditions.

    The ZeroD model uses the underlying surface kinetics (chem) to compute the
    saturation curve, which describes how surface coverage evolves with time
    under constant exposure conditions.

    Parameters
    ----------
    chem : SurfaceKinetics
        The surface kinetics model (e.g., ALDideal, ALDsoft, ALDchem) that defines
        the sticking probabilities, recombination, and other reaction parameters.
    T : float
        Temperature in Kelvin.
    p : float
        Precursor partial pressure in Pascals.

    Attributes
    ----------
    chem : SurfaceKinetics
        The surface kinetics model.
    T : float
        Operating temperature in Kelvin.
    p : float
        Precursor partial pressure in Pascals.
    vth : float
        Thermal velocity of the precursor at the current temperature (inherited from DoseModel).

    Methods
    -------
    saturation_curve(T=None, p=None)
        Generate the time-dependent saturation curve at the specified conditions.

    Examples
    --------
    >>> from aldsim.chem import Precursor, ALDideal
    >>> tma = Precursor('TMA', mass=144.17)
    >>> kinetics = ALDideal(prec=tma, nsites=1e19, p_stick=0.1)
    >>> model = ZeroD(chem=kinetics, T=473.15, p=100)
    >>> time, coverage = model.saturation_curve()
    """

    def __init__(self, chem, T, p):
        super().__init__(chem, T, p)

    def saturation_curve(self, T=None, p=None):
        """
        Generate the time-dependent saturation curve.

        Computes the coverage evolution over time for the self-limited process
        at the specified (or current) temperature and pressure conditions.

        Parameters
        ----------
        T : float, optional
            Temperature in Kelvin. If provided, updates the model's temperature.
            If None, uses the current temperature.
        p : float, optional
            Precursor partial pressure in Pascals. If provided, updates the
            model's pressure. If None, uses the current pressure.

        Returns
        -------
        tuple of ndarray
            (time, coverage) arrays where:
            - time: Array of time points in seconds
            - coverage: Array of corresponding coverage values (0 to 1)
        """
        if T is not None:
            self.T = T
        if p is not None:
            self.p = p
        return self.chem.saturation_curve(self.T, self.p)

