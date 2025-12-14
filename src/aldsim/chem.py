#Copyright © 2024-2025, UChicago Argonne, LLC

from .aldutils import calc_vth
from .constants import kb, Rgas, Nav

import numpy as np

_precursor_mass = {
    'TMA' : 144.17,
    'H2O' : 18.01,
    'O': 16,
    'N': 14
}

class Precursor:
    """
    Defines a precursor molecule
    """

    def __init__(self, name='None', mass=None, ligands=None):
        self.name = name
        if mass is None:
            self.mass = _precursor_mass[self.name]
        else:
            self.mass = mass
        self.ligands = ligands

    def vth(self, T):
        """Calculate the mean thermal velocity at temperature T (in K)"""
        return calc_vth(self.mass, T)
    
    def Jwall(self, T, p, in_mols=False):
        """Calculate the flux per unit area for a given temperature (in K) and pressure (in Pa)
        
        If in_mols is True returns the mols per surface area, otherwise it returns the molecules per
        surface area impinging on a surface
        """
        if in_mols:
            return 0.25*self.vth(T)*p/(Rgas*T)
        else:
            return 0.25*self.vth(T)*p/(kb*T)


class SurfaceKinetics:
    """
    Base class for self-limited kinetics

    It assumes that a fraction of the surface is reactive and comprised of reaction
    sites of equal surface area. The default is that all the surface is reactive.

    """

    def __init__(self, prec, nsites, f=1):
        self.prec = prec
        self._f = f
        self.nsites = nsites

    @property
    def site_area(self):
        """Area of a single reaction site"""
        return self._s0
    
    @site_area.setter
    def site_area(self, value):
        self._s0 = value
        self._nsites = self._f/self._s0
    
    @property
    def nsites(self):
        """Number of reactive sites per surface area"""
        return self._nsites

    @nsites.setter
    def nsites(self, value):
        self._nsites = value
        self._s0 = self._f/self._nsites

    @property
    def nsites_mol(self):
        """Number of reactive sites per surface area in mols"""
        return self.nsites/Nav
    
    @property
    def f(self):
        """Fraction of reactive sites"""
        return self._f 

    @f.setter
    def f(self, value):
        self._f = value
        self._nsites = self._f/self._s0

    def sticking_prob(self, *args):
        pass

    def sticking_prob_av(self, *args):
        pass

    def vth(self, T):
        return self.prec.vth(T)
    
    def Jwall(self, T, p):
        return self.prec.Jwall(T, p)


class ALDideal(SurfaceKinetics):
    """
    Ideal first-order irreversible Langmuir kinetics for atomic layer deposition.

    This class implements a surface kinetics model for ALD processes where precursor
    molecules undergo chemisorption on reactive surface sites. The model accounts for
    both direct sticking and recombination pathways, with coverage-dependent probabilities.

    The model assumes:
    - First-order Langmuir adsorption kinetics
    - Self-limiting surface reactions
    - Uniform reactive sites with equal surface area
    - Coverage-dependent recombination, with different values depending on the state of the site

    Parameters
    ----------
    prec : Precursor
        The precursor molecule for this half-reaction. Contains properties like
        molecular mass needed to calculate thermal velocity and wall flux.
    nsites : float
        Number of reactive sites per unit surface area (sites/m²). This determines
        the site area via the relationship: site_area = f / nsites.
    p_stick : float
        Sticking probability (0 ≤ p_stick ≤ 1) for precursor molecules
        interacting with bare reactive sites. Represents the probability that a
        gas-phase molecule chemisorbs upon collision with an unreacted site.
    p_rec0 : float, optional
        Recombination probability on bare sites (default: 0). This is the probability
        that a precursor molecule reacts via recombination with the bare surface, for
        instance in the case of plasma species.
    p_rec1 : float, optional
        Recombination probability on fully reacted sites (default: 0). This is the
        probability of recombination when the surface is saturated (coverage = 1).
        The actual recombination probability varies linearly with coverage between
        p_rec0 and p_rec1.
    f : float, optional
        Fraction of the surface that is reactive (0 < f ≤ 1, default: 1). For f < 1,
        only a portion of the surface contains reactive sites.
    dm : float, optional
        mass gain during the cycle (default: 1). Can be used to generate qcm data.

    Attributes
    ----------
    name : str
        Identifier for this kinetics model type ('ideal')
    p_stick0 : float
        Initial sticking probability
    p_rec0 : float
        Recombination probability on bare sites
    p_rec1 : float
        Recombination probability on reacted sites
    dm : float
        mass gain

    Methods
    -------
    sticking_prob(cov=0)
        Calculate the coverage-dependent sticking probability
    recomb_prob(cov=0)
        Calculate the coverage-dependent recombination probability
    react_prob(cov=0)
        Calculate the total reaction probability (sticking + recombination)
    t0(T, p)
        Calculate the characteristic saturation time
    saturation_curve(T, p)
        Generate the time-dependent saturation curve

    Examples
    --------
    >>> tma = Precursor('TMA', mass=144.17)
    >>> kinetics = ALDideal(prec=tma, nsites=1e19, p_stick=0.1, p_rec0=0.01, p_rec1=0.05)
    >>> # Calculate sticking probability at 50% coverage
    >>> s = kinetics.sticking_prob(cov=0.5)
    >>> # Get characteristic time at 200°C and 100 Pa
    >>> t_sat = kinetics.t0(T=473.15, p=100)
    """

    name = 'ideal'

    def __init__(self, prec, nsites, p_stick, p_rec0=0, p_rec1=0, f=1, dm=1):
        self.p_stick0 = p_stick
        self.p_rec0 = p_rec0
        self.p_rec1 = p_rec1
        self.dm = dm
        super().__init__(prec, nsites, f)

    def sticking_prob(self, cov=0):
        """
        Calculate the coverage-dependent sticking probability.

        The sticking probability decreases linearly with coverage as reactive sites
        become occupied: S(θ) = f * p_stick0 * (1 - θ)

        Parameters
        ----------
        cov : float, optional
            Surface coverage (0 ≤ cov ≤ 1, default: 0). Represents the fraction of
            reactive sites that have already reacted with precursor molecules.

        Returns
        -------
        float
            Sticking probability at the given coverage (0 ≤ S ≤ f*p_stick0)
        """
        return self.f*self.p_stick0*(1-cov)

    def recomb_prob(self, cov=0):
        """
        Calculate the coverage-dependent recombination probability.

        The recombination probability varies linearly between bare and reacted sites:
        R(θ) = p_rec0 + f * θ * (p_rec1 - p_rec0)

        This accounts for ligand recombination reactions that release volatile species
        without necessarily contributing to film growth.

        Parameters
        ----------
        cov : float, optional
            Surface coverage (0 ≤ cov ≤ 1, default: 0)

        Returns
        -------
        float
            Recombination probability at the given coverage
        """
        return self.p_rec0 + self.f*cov*(self.p_rec1-self.p_rec0)

    def react_prob(self, cov=0):
        """
        Calculate the total reaction probability.

        The total reaction probability is the sum of sticking and recombination
        probabilities: P_react(θ) = S(θ) + R(θ)

        This represents the overall probability that a gas-phase precursor molecule
        undergoes any surface reaction upon collision.

        Parameters
        ----------
        cov : float, optional
            Surface coverage (0 ≤ cov ≤ 1, default: 0)

        Returns
        -------
        float
            Total reaction probability at the given coverage
        """
        return self.sticking_prob(cov) + self.recomb_prob(cov)

    def sticking_prob_av(self, av):
        """
        Calculate the average sticking probability based on the fraction of available sites.

        This method computes the sticking probability based on site availability
        rather than coverage. Useful for certain averaging schemes in simulations.

        Parameters
        ----------
        av : float
            Fraction of available sites (0 ≤ av ≤ 1). Represents the fraction of sites
            available for reaction.

        Returns
        -------
        float
            Average sticking probability
        """
        return self.f*self.p_stick0*av

    def t0(self, T, p):
        """
        Calculate the characteristic saturation time.

        The characteristic time represents the timescale for surface saturation
        under ideal first-order Langmuir kinetics:
        t0 = 1 / (site_area * J_wall * p_stick0)

        where J_wall is the molecular flux to the wall.

        Parameters
        ----------
        T : float
            Temperature in Kelvin
        p : float
            Precursor partial pressure in Pascals

        Returns
        -------
        float
            Characteristic saturation time in seconds
        """
        return 1.0/(self.site_area*self.Jwall(T, p)*self.p_stick0)

    def saturation_curve(self, T, p):
        """
        Generate the time-dependent saturation curve.

        Calculates the coverage evolution over time for ideal first-order Langmuir
        kinetics: θ(t) = 1 - exp(-t/t0)

        The method automatically determines an appropriate time range based on
        the characteristic saturation time.

        Parameters
        ----------
        T : float
            Temperature in Kelvin
        p : float
            Precursor partial pressure in Pascals

        Returns
        -------
        tuple of ndarray
            (time, coverage) arrays where:
            - time: Array of time points in seconds
            - coverage: Array of corresponding coverage values (0 to 1)
        """
        t0 = self.t0(T,p)
        tscale = 5*t0
        logtscale = np.log10(tscale)
        scale = int(logtscale)
        if logtscale < 0:
            scale -= 1
        factor = int(10**(logtscale-scale))+1
        tmax= factor*10**scale
        print(t0, tmax)
        dt = tmax/100
        x = np.arange(0, tmax, dt)
        y = 1-np.exp(-x/t0)
        return x, y


class ALDsoft(SurfaceKinetics):
    """
    First-order irreversible Langmuir kinetics with two independent reaction pathways.

    This class extends the ideal ALD kinetics model to account for surface heterogeneity
    by considering two distinct types of reactive sites, each with different sticking
    probabilities. This dual-pathway model is useful for describing "soft saturation"
    behavior where different surface sites saturate at different rates.

    The model assumes:
    - Two independent populations of reactive sites (pathway 1 and pathway 2)
    - First-order Langmuir adsorption kinetics for each pathway
    - Self-limiting surface reactions on both pathways
    - Each pathway has uniform reactive sites with equal surface area
    - The fraction of unreactive sites is given by: f_unreactive = 1 - (f1 + f2)

    The main difference from ALDideal is that ALDsoft considers two separate reaction
    pathways with independent coverage evolution. f1 and f2 are the fractions of the
    surface sites for each pathway. Normally f1 + f2 = 1 unless there is a fraction
    of unreactive sites.

    Parameters
    ----------
    prec : Precursor
        The precursor molecule for this half-reaction. Contains properties like
        molecular mass needed to calculate thermal velocity and wall flux.
    nsites : float
        Number of reactive sites per unit surface area (sites/m²) for each pathway.
        This determines the site area via the relationship: site_area = (f1 + f2) / nsites.
    p_stick1 : float
        Sticking probability (0 ≤ p_stick1 ≤ 1) for precursor molecules interacting
        with bare reactive sites on pathway 1. Represents the probability that a
        gas-phase molecule chemisorbs upon collision with an unreacted site of type 1.
    p_stick2 : float
        Sticking probability (0 ≤ p_stick2 ≤ 1) for precursor molecules interacting
        with bare reactive sites on pathway 2. Typically p_stick2 ≠ p_stick1,
        leading to different saturation timescales for each pathway.
    f1 : float
        Fraction of the surface covered by pathway 1 reactive sites (0 < f1 ≤ 1).
        This represents the population of sites with sticking probability p_stick1.
    f2 : float, optional
        Fraction of the surface covered by pathway 2 reactive sites (0 < f2 ≤ 1).
        If not provided, defaults to (1 - f1), assuming no unreactive sites.
        For surfaces with unreactive sites, set f2 explicitly such that f1 + f2 < 1.

    Attributes
    ----------
    name : str
        Identifier for this kinetics model type ('softsat')
    p_stick1 : float
        Sticking probability for pathway 1 sites
    p_stick2 : float
        Sticking probability for pathway 2 sites
    f1 : float
        Fraction of pathway 1 sites
    f2 : float
        Fraction of pathway 2 sites

    Methods
    -------
    sticking_prob(cov1=0, cov2=0)
        Calculate the total coverage-dependent sticking probability from both pathways
    sticking_prob_av(av1, av2)
        Calculate the average sticking probability based on site availability
    t0(T, p)
        Calculate the characteristic saturation times for both pathways
    saturation_curve(T, p)
        Generate the time-dependent saturation curve combining both pathways

    Examples
    --------
    >>> tma = Precursor('TMA', mass=144.17)
    >>> # Surface with 60% fast-reacting sites and 40% slow-reacting sites
    >>> kinetics = ALDsoft(prec=tma, nsites=1e19, p_stick1=0.2, p_stick2=0.02,
    ...                    f1=0.6, f2=0.4)
    >>> # Calculate sticking probability with 50% coverage on pathway 1, 20% on pathway 2
    >>> s = kinetics.sticking_prob(cov1=0.5, cov2=0.2)
    >>> # Get characteristic times for both pathways at 200°C and 100 Pa
    >>> t1, t2 = kinetics.t0(T=473.15, p=100)
    >>>
    >>> # Example with unreactive sites (10% of surface is unreactive)
    >>> kinetics_partial = ALDsoft(prec=tma, nsites=1e19, p_stick1=0.2,
    ...                            p_stick2=0.05, f1=0.7, f2=0.2)
    >>> # Here f1 + f2 = 0.9, meaning 10% of the surface is unreactive
    """

    name = 'softsat'

    def __init__(self, prec, nsites, p_stick1, p_stick2, f1, f2=None):
        self.p_stick1 = p_stick1
        self.p_stick2 = p_stick2
        self.f1 = f1
        if f2 is None:
            self.f2 = 1-self.f1
        else:
            self.f2 = f2
        super().__init__(prec, nsites, self.f1+self.f2)

    def sticking_prob(self, cov1=0, cov2=0):
        """
        Calculate the total coverage-dependent sticking probability from both pathways.

        The total sticking probability is the sum of contributions from both pathways,
        each decreasing linearly with their respective coverage:
        S_total(θ1, θ2) = f1 * p_stick1 * (1 - θ1) + f2 * p_stick2 * (1 - θ2)

        This allows the two pathways to saturate independently at different rates,
        producing the characteristic "soft saturation" behavior.

        Parameters
        ----------
        cov1 : float, optional
            Surface coverage on pathway 1 sites (0 ≤ cov1 ≤ 1, default: 0).
            Represents the fraction of pathway 1 reactive sites that have reacted.
        cov2 : float, optional
            Surface coverage on pathway 2 sites (0 ≤ cov2 ≤ 1, default: 0).
            Represents the fraction of pathway 2 reactive sites that have reacted.

        Returns
        -------
        float
            Total sticking probability at the given coverages
        """
        return self.f1*self.p_stick1*(1-cov1) + self.f2*self.p_stick2*(1-cov2)

    def sticking_prob_av(self, av1, av2):
        """
        Calculate the average sticking probability based on site availability for both pathways.

        This method computes the total sticking probability based on site availability
        rather than coverage for each pathway. Useful for certain averaging schemes
        in simulations.

        Parameters
        ----------
        av1 : float
            Fraction of available sites on pathway 1 (0 ≤ av1 ≤ 1). Represents the
            fraction of pathway 1 sites available for reaction.
        av2 : float
            Fraction of available sites on pathway 2 (0 ≤ av2 ≤ 1). Represents the
            fraction of pathway 2 sites available for reaction.

        Returns
        -------
        float
            Total average sticking probability from both pathways
        """
        return self.f1*self.p_stick1*av1 + self.f2*self.p_stick2*av2

    def t0(self, T, p):
        """
        Calculate the characteristic saturation times for both pathways.

        Each pathway has its own characteristic time representing the timescale for
        saturation under ideal first-order Langmuir kinetics:
        t0_i = 1 / (site_area * J_wall * p_stick_i)

        where J_wall is the molecular flux to the wall and i = 1 or 2.

        The difference between t1 and t2 determines the extent of soft saturation
        behavior. Large differences lead to more pronounced multi-step saturation curves.

        Parameters
        ----------
        T : float
            Temperature in Kelvin
        p : float
            Precursor partial pressure in Pascals

        Returns
        -------
        tuple of float
            (t1, t2) where:
            - t1: Characteristic saturation time for pathway 1 in seconds
            - t2: Characteristic saturation time for pathway 2 in seconds
        """
        t1 = 1.0/(self.site_area*self.Jwall(T, p)*self.p_stick1)
        t2 = 1.0/(self.site_area*self.Jwall(T, p)*self.p_stick2)
        return t1, t2

    def saturation_curve(self, T, p):
        """
        Generate the time-dependent saturation curve combining both pathways.

        Calculates the total coverage evolution over time as a weighted sum of two
        independent first-order Langmuir processes:
        θ_total(t) = [f1*(1 - exp(-t/t1)) + f2*(1 - exp(-t/t2))] / (f1 + f2)

        The resulting curve exhibits soft saturation behavior when t1 ≠ t2, with an
        initial fast saturation phase (dominated by the faster pathway) followed by
        a slower approach to complete saturation.

        The method automatically determines an appropriate time range based on
        the slower of the two characteristic saturation times.

        Parameters
        ----------
        T : float
            Temperature in Kelvin
        p : float
            Precursor partial pressure in Pascals

        Returns
        -------
        tuple of ndarray
            (time, coverage) arrays where:
            - time: Array of time points in seconds
            - coverage: Array of corresponding total coverage values (0 to 1)
        """
        t1, t2 = self.t0(T,p)
        t0 = max(t1, t2)
        tscale = 5*t0
        logtscale = np.log10(tscale)
        scale = int(logtscale)
        if logtscale < 0:
            scale -= 1
        factor = int(10**(logtscale-scale))+1
        tmax= factor*10**scale
        print(t0, tmax)
        dt = tmax/100
        x = np.arange(0, tmax, dt)
        y = (self.f1*(1-np.exp(-x/t1)) + self.f2*(1-np.exp(-x/t2)))/(self.f1+self.f2)
        return x, y



class ALDProcess:

    def __init__(self, chem1, chem2, n1, n2, T=None):

        self.chem1 = chem1
        self.n1 = n1
        self.n2 = n2
        self.chem2 = chem2
        self.T = T

    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, value):
        self._T = value
        self.chem1.T = value
        self.chem2.T = value

    
    


    





