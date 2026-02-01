#Copyright © 2024-Present, UChicago Argonne, LLC

"""ALD transport in high aspect ratio features"""

from .base import DoseModel
from aldsim.constants import kb
from aldsim.core.diffusion import DiffusionViaND

import numpy as np

class KnudsenVia(DoseModel):
    """
    Model of ALD particle coating in high aspect ratio circular via
    ----------

    chem : ALDchem
        Surface chemistry object defining the reaction kinetics. Must have
        a single reaction pathway (single_path=True).
    p : float
        Precursor partial pressure in Pa.
    T : float
        Temperature in K.
    AR : float
        Aspect ratio of the via

    """

    def __init__(self, chem, p, T, AR):
        if not chem.single_path:
            raise NotImplementedError("KnudsenVia only supports single pathway chemistry")
        super().__init__(chem, T, p)
        self.AR = AR
        self.base_model = DiffusionViaND(p_stick0=self.chem.p_stick1,
            p_rec0=self.chem.p_rec0, p_rec1=self.chem.p_rec1, AR=self.AR)

    def t0(self):
        """Calculate the characteristic time scale for the ALD process.

        Returns
        -------
        float
            Characteristic time in seconds, based on the precursor flux,
            sticking probability, and surface site density.
        """
        return 4*kb*self.T/(self.chem.site_area*self.vth*self.p*self.chem.p_stick1)

    @property
    def z(self):
        """Spatial coordinates along the via depth.

        Returns
        -------
        ndarray
            Array of z positions (normalized to the via diameter) at the
            center of each discretized segment, from the via entrance to
            the bottom.
        """
        return np.array([(0.5+i)*self.base_model.dz
            for i in range(self.AR*self.base_model.nsegments)])

    def saturation_curve(self):
        """Run simulation until full surface saturation is reached.

        Executes the diffusion-reaction model until the mean coverage
        reaches 99%, saving coverage profiles at regular intervals.

        Returns
        -------
        times : list of float
            List of times in seconds corresponding to saved coverage profiles.
        z : ndarray
            Spatial coordinates along the via depth, normalized to the
            via diameter.
        coverage : list of ndarray
            List of coverage arrays at saved time points. Each array
            represents the coverage profile along the via depth.
        """
        t, cov = self.base_model.run_until_cov()
        t0 = self.t0()
        tout = [ti*t0 for ti in t]
        return tout, self.z, cov
    
    def run(self, tdose=None, save_every=None, dtrun=None):
        """Run ALD simulation in the via.

        Parameters
        ----------
        tdose : float, optional
            Dose time in seconds. If None (default), runs until saturation
            and returns the same result as saturation_curve().
        save_every : float, optional
            Time interval at which coverage profiles are saved (in normalized
            units). If None, defaults to tdose/5.
        dtrun : float, optional
            Time step for numerical integration (in normalized units).
            If None, uses the default from the base model.

        Returns
        -------
        times : list of float
            List of times in seconds corresponding to saved coverage profiles.
        z : ndarray
            Spatial coordinates along the via depth, normalized to the
            via diameter.
        coverage : list of ndarray
            List of coverage arrays at saved time points. Each array
            represents the coverage profile along the via depth.
        """
        if tdose is None:
            return self.saturation_curve()
        else:
            trun = tdose/self.t0()
            if save_every is None:
                save_every = trun/5
            if dtrun is None:
                tlist, covlist = self.base_model.run(max_time=trun, save_every=save_every)
            else:
                tlist, covlist = self.base_model.run(max_time=trun, save_every=save_every, dt=dtrun)

        t0 = self.t0()
        tout = [ti*t0 for ti in tlist]
        return tout, self.z, covlist
