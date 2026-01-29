#Copyright © 2024-Present, UChicago Argonne, LLC

import pytest
from aldsim.models.batch import FluidizedBed
from aldsim.chem import Precursor, ALDchem


class TestFluidizedBed:

    def setup_method(self):
        """Set up test fixtures."""
        self.prec = Precursor(mass=100)
        self.chem = ALDchem(self.prec, nsites=1e19, p_stick1=0.01)
        # Typical test parameters
        self.p = 100       # Pa, precursor partial pressure
        self.p0 = 1e5      # Pa, carrier gas pressure
        self.T = 450       # K, temperature
        self.S = 1.0       # m², total surface area
        self.flow = 100    # sccm, gas flow rate

    def test_init(self):
        """Test basic initialization."""
        drum = FluidizedBed(self.chem, self.p, self.p0, self.T, self.S, self.flow)
        assert drum.S == self.S
        assert drum.p0 == self.p0
        assert drum.flow0 == self.flow
        assert drum.T == self.T
        assert drum.p == self.p

    def test_init_rejects_dual_pathway(self):
        """Test that dual pathway chemistry raises NotImplementedError."""
        dual_chem = ALDchem(self.prec, nsites=1e19, p_stick1=0.01, p_stick2=0.001, f2=0.2)
        with pytest.raises(NotImplementedError):
            FluidizedBed(dual_chem, self.p, self.p0, self.T, self.S, self.flow)

    def test_flow_positive(self):
        """Test that flow() returns a positive value."""
        drum = FluidizedBed(self.chem, self.p, self.p0, self.T, self.S, self.flow)
        assert drum.flow() > 0

    def test_da_positive(self):
        """Test that Da() returns a positive Damkohler number."""
        drum = FluidizedBed(self.chem, self.p, self.p0, self.T, self.S, self.flow)
        da = drum.Da()
        assert da > 0

    def test_t0_positive(self):
        """Test that t0() returns a positive characteristic time."""
        drum = FluidizedBed(self.chem, self.p, self.p0, self.T, self.S, self.flow)
        t0 = drum.t0()
        assert t0 > 0

    def test_saturation_curve(self):
        """Test saturation_curve returns proper arrays."""
        drum = FluidizedBed(self.chem, self.p, self.p0, self.T, self.S, self.flow)
        t, cov = drum.saturation_curve()
        # Arrays should have same shape
        assert t.shape == cov.shape
        # Time should be monotonically increasing
        assert all(t[i] < t[i+1] for i in range(len(t)-1))
        # Coverage should start near zero and end near one
        assert cov[0] < 0.1
        assert cov[-1] > 0.9

    def test_run_default(self):
        """Test run() with default parameters."""
        drum = FluidizedBed(self.chem, self.p, self.p0, self.T, self.S, self.flow)
        t, cov, x = drum.run()
        # Arrays should have same shape
        assert t.shape == cov.shape
        assert t.shape == x.shape
        # Coverage should start near zero and end near one
        assert cov[0] < 0.1
        assert cov[-1] > 0.9

    def test_run_with_tdose(self):
        """Test run() with specified dose time."""
        drum = FluidizedBed(self.chem, self.p, self.p0, self.T, self.S, self.flow)
        tdose = drum.t0() * 2  # Run for 2x characteristic time
        t, cov, x = drum.run(tdose=tdose)
        assert t.shape == cov.shape
        assert t.shape == x.shape

    def test_run_with_tdose_and_dt(self):
        """Test run() with specified dose time and time step."""
        drum = FluidizedBed(self.chem, self.p, self.p0, self.T, self.S, self.flow)
        tdose = drum.t0() * 2
        dt = drum.t0() * 0.1
        t, cov, x = drum.run(tdose=tdose, dt=dt)
        assert t.shape == cov.shape
        assert t.shape == x.shape
