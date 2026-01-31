#Copyright © 2024-Present, UChicago Argonne, LLC

import pytest
from aldsim.models.spatial import SpatialWellMixed
from aldsim.chem import Precursor, ALDchem


class TestSpatialWellMixed:

    def setup_method(self):
        """Set up test fixtures."""
        self.prec = Precursor(mass=100)
        self.chem = ALDchem(self.prec, nsites=1e19, p_stick1=0.01)
        # Typical test parameters
        self.p = 100       # Pa, precursor partial pressure
        self.p0 = 1e5      # Pa, carrier gas pressure
        self.T = 450       # K, temperature
        self.L = 0.02      # m, length of spatial ALD zone
        self.W = 0.1       # m, width of spatial ALD zone
        self.S = 1.0       # m², total surface area
        self.flow = 100    # sccm, gas flow rate

    def test_init(self):
        """Test basic initialization."""
        model = SpatialWellMixed(self.chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)
        assert model.S == self.S
        assert model.p0 == self.p0
        assert model.flow0 == self.flow
        assert model.T == self.T
        assert model.p == self.p

    def test_init_rejects_dual_pathway(self):
        """Test that dual pathway chemistry raises NotImplementedError."""
        dual_chem = ALDchem(self.prec, nsites=1e19, p_stick1=0.01, p_stick2=0.001, f2=0.2)
        with pytest.raises(NotImplementedError):
            SpatialWellMixed(dual_chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)

    def test_flow_positive(self):
        """Test that flow() returns a positive value."""
        model = SpatialWellMixed(self.chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)
        assert model.flow() > 0

    def test_da_positive(self):
        """Test that Da() returns a positive Damkohler number."""
        model = SpatialWellMixed(self.chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)
        da = model.Da()
        assert da > 0

    def test_t0_positive(self):
        """Test that t0() returns a positive characteristic time."""
        model = SpatialWellMixed(self.chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)
        t0 = model.t0()
        assert t0 > 0

    def test_saturation_curve(self):
        """Test saturation_curve returns proper arrays."""
        model = SpatialWellMixed(self.chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)
        velocity, cov = model.saturation_curve()
        # Arrays should have same shape
        assert velocity.shape == cov.shape
        # Velocity should be monotonically decreasing (high speed = low coverage)
        assert all(velocity[i] > velocity[i+1] for i in range(len(velocity)-1))
        # Coverage should start near zero (high velocity) and end near one (low velocity)
        assert cov[0] < 0.1
        assert cov[-1] > 0.9

    def test_run_default(self):
        """Test run() with default parameters."""
        model = SpatialWellMixed(self.chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)
        velocity, cov, x = model.run()
        # Arrays should have same shape
        assert velocity.shape == cov.shape
        assert velocity.shape == x.shape
        # Coverage should start near zero (high velocity) and end near one (low velocity)
        assert cov[0] < 0.1
        assert cov[-1] > 0.9

    def test_run_with_umax(self):
        """Test run() with specified maximum velocity."""
        model = SpatialWellMixed(self.chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)
        umax = 0.1  # m/s maximum web velocity
        velocity, cov, x = model.run(umax=umax)
        assert velocity.shape == cov.shape
        assert velocity.shape == x.shape

    def test_run_with_umax_and_du(self):
        """Test run() with specified maximum velocity and velocity step."""
        model = SpatialWellMixed(self.chem, self.p, self.p0, self.T, self.flow, self.L, self.W, self.S)
        umax = 0.1  # m/s maximum web velocity
        du = 0.001  # m/s velocity step
        velocity, cov, x = model.run(umax=umax, du=du)
        assert velocity.shape == cov.shape
        assert velocity.shape == x.shape
