#Copyright © 2024-Present, UChicago Argonne, LLC

"""Batch particle coating ALD models"""

from .base import DoseModel
from aldsim.constants import kb
from aldsim.core.particle import WellMixedParticleND
from aldsim.core.particle.fluidizedbed import FluidizedBedND


class RotatingDrum(DoseModel):
    """
    Model of ALD particle coating in a rotating drum or well-stirred reactor.

    This model simulates the coating of particles in a batch reactor where
    particles are well-mixed, such as a rotating drum reactor. It assumes
    uniform precursor concentration throughout the reactor volume and
    first-order Langmuir kinetics for the surface reactions.

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
    S : float
        Total surface area of the particles to be coated in m².
    flow : float
        Gas flow rate to the reactor in sccm (standard cubic centimeters per minute).

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
        Run the simulation and return time, coverage, and precursor concentration.

    Raises
    ------
    NotImplementedError
        If chem has more than one reaction pathway.
    """

    def __init__(self, chem, p, p0, T, S, flow):
        if not chem.single_path:
            raise NotImplementedError("RotatingDrum only supports single pathway chemistry")
        super().__init__(chem, T, p)
        self.S = S
        self.p0 = p0
        self.flow0 = flow
        da = self.Da()
        self.base_model = WellMixedParticleND(da)

    def flow(self):
        return (1e-6*self.flow0/60)*1e5/self.p0*(self.T/300)

    def Da(self):
        flow = self.flow()
        return 0.25*self.S/flow*self.chem.sticking_prob()*self.vth

    def t0(self):
        return kb*self.T*self.S/(self.flow()*self.chem.site_area*self.p)

    def saturation_curve(self):
        self.base_model.Da = self.Da()
        t, cov = self.base_model.saturation_curve()
        return t*self.t0(), cov
    
    def run(self, tdose=None, dt=None):
        self.base_model.Da = self.Da()
        if tdose is None:
            t, cov, x = self.base_model.run()

        else:
            trun = tdose/self.t0()
            if dt is None:
                dtrun = 0.01
            else:
                dtrun = dt/self.t0()
            t, cov, x = self.base_model.run(trun, dtrun)
            
        return t*self.t0(), cov, x



class FluidizedBed(DoseModel):
    """
    Model of ALD particle coating in a fluidized bed reactor.

    This model simulates the coating of particles in a batch fluidized bed reactor
    where particles are well-mixed and precursor transport follows plug flow
    approximation. It assumes first-order Langmuir kinetics for the surface reactions.

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
    S : float
        Total surface area of the particles to be coated in m².
    flow : float
        Gas flow rate to the reactor in sccm (standard cubic centimeters per minute).

    Attributes
    ----------
    S : float
        Total particle surface area (m²).
    p0 : float
        Carrier gas pressure (Pa).
    flow0 : float
        Gas flow rate (sccm).
    base_model : FluidizedBedND
        Underlying dimensionless fluidized bed model.

    Methods
    -------
    Da()
        Calculate the Damkohler number for the reactor.
    t0()
        Calculate the characteristic saturation time.
    saturation_curve()
        Generate the time-dependent saturation curve.
    run(tdose=None, dt=None)
        Run the simulation and return time, coverage, and precursor concentration.

    Raises
    ------
    NotImplementedError
        If chem has more than one reaction pathway.
    """

    def __init__(self, chem, p, p0, T, S, flow):
        if not chem.single_path:
            raise NotImplementedError("FluidizedBed only supports single pathway chemistry")
        super().__init__(chem, T, p)
        self.S = S
        self.p0 = p0
        self.flow0 = flow
        da = self.Da()
        self.base_model = FluidizedBedND(da)

    def flow(self):
        return (1e-6*self.flow0/60)*1e5/self.p0*(self.T/300)

    def Da(self):
        flow = self.flow()
        return 0.25*self.S/flow*self.chem.sticking_prob()*self.vth

    def t0(self):
        return kb*self.T*self.S/(self.flow()*self.chem.site_area*self.p)

    def saturation_curve(self):
        self.base_model.Da = self.Da()
        t, cov = self.base_model.saturation_curve()
        return t*self.t0(), cov

    def run(self, tdose=None, dt=None):
        self.base_model.Da = self.Da()
        if tdose is None:
            t, cov, x = self.base_model.run()

        else:
            trun = tdose/self.t0()
            if dt is None:
                dtrun = 0.01
            else:
                dtrun = dt/self.t0()
            t, cov, x = self.base_model.run(trun, dtrun)

        return t*self.t0(), cov, x

