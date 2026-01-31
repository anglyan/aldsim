from .base import DoseModel
from aldsim.constants import kb
from aldsim.core.particle import SpatialWellMixedND


class SpatialWellMixed(DoseModel):
    """
    Model of ALD coating in a spatial ALD reactor using a well-stirred precursor approximation.

    This model simulates ALD coating in a spatial reactor where the precursor
    is well-mixed in the direction perpendicular to the substrate movement.
    It assumes uniform precursor concentration throughout the reactor volume
    and first-order Langmuir kinetics for the surface reactions.

    The model supports two configurations:
    - Flat surface coating: When S is None, the model treats the substrate as
      a flat surface with area L × W.
    - Particle bed coating: When S is provided, it represents the total surface
      area of particles to be coated within the reactor zone.

    Parameters
    ----------
    chem : ALDchem
        Surface chemistry object defining the reaction kinetics. Must have
        a single reaction pathway (single_path=True).
    p : float
        Precursor partial pressure in Pa.
    p0 : float
        Carrier gas pressure in Pa.
    T : float
        Temperature in K.
    flow : float
        Gas flow rate to the reactor in sccm (standard cubic centimeters per minute).
    L : float
        Length of the spatial ALD zone (in m).
    W : float
        Width of the spatial ALD zone (in m).
    S : float, optional
        Total surface area to be coated in m². If None (default), the surface
        area is calculated as L × W for a flat substrate.

    Attributes
    ----------
    S : float
        Total particle surface area (m²).
    p0 : float
        Carrier gas pressure (Pa).
    flow0 : float
        Gas flow rate (sccm).
    base_model : WellMixedParticleND
        Underlying dimensionless well-mixed particle model.

    Methods
    -------
    Da()
        Calculate the Damkohler number for the reactor.
    t0()
        Calculate the characteristic saturation time.
    saturation_curve()
        Generate the time-dependent saturation curve.
    run(tdose=None, dt=None)
        Run the simulation and return residence time, coverage, and precursor concentration.

    Raises
    ------
    NotImplementedError
        If chem has more than one reaction pathway.
    """

    def __init__(self, chem, p, p0, T, flow, L, W, S=None):
        if not chem.single_path:
            raise NotImplementedError("SpatialWellMixed only supports single pathway chemistry")
        super().__init__(chem, T, p)
        self.W = W
        self.L = L
        if S is None:
            self.S = W*L
        else:
            self.S = S
        self.p0 = p0
        self.flow0 = flow
        da = self.Da()
        self.base_model = SpatialWellMixedND(da)

    def flow(self):
        return (1e-6*self.flow0/60)*1e5/self.p0*(self.T/300)

    def Da(self):
        flow = self.flow()
        return 0.25*self.S/flow*self.chem.sticking_prob()*self.vth

    def t0(self):
        return kb*self.T*self.S/(self.flow()*self.chem.site_area*self.p)

    def saturation_curve(self):
        """Return the saturation curve as a function of the web velocity"""
        
        self.base_model.Da = self.Da()
        tmax = max(5, 10/self.base_model.Da)
        t, cov = self.base_model.saturation_curve(tmax=tmax)

        if t[0] == 0:
            t = t[1:]
            cov = cov[1:]
        return self.L/(t*self.t0()), cov
    
    def run(self, umax=None, du=None):
        self.base_model.Da = self.Da()
        if umax is None:
            t, cov, x = self.base_model.run()

        else:
            if du is None:
                du = 0.001*umax
            trun = self.L/du/self.t0()
            dtrun = self.L/(umax*self.t0())
            t, cov, x = self.base_model.run(trun, dtrun)
            
        return self.L/(t*self.t0()), cov, x

