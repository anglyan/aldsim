#Copyright © 2024-Present, UChicago Argonne, LLC

import pytest
from aldsim.models.knudsen import KnudsenVia
from aldsim.chem import Precursor, ALDchem


class TestKnudsenVia:

    def setup_method(self):
        """Set up test fixtures."""
        self.prec = Precursor(mass=100)
        self.chem = ALDchem(self.prec, nsites=1e19, p_stick1=0.01)
        # Typical test parameters
        self.p = 100       # Pa, precursor partial pressure
        self.T = 450       # K, temperature
        self.AR = 10       # Aspect ratio

    def test_init(self):
        """Test basic initialization."""
        model = KnudsenVia(self.chem, self.p, self.T, self.AR)
        assert model.AR == self.AR
        assert model.T == self.T
        assert model.p == self.p

    def test_init_rejects_dual_pathway(self):
        """Test that dual pathway chemistry raises NotImplementedError."""
        dual_chem = ALDchem(self.prec, nsites=1e19, p_stick1=0.01, p_stick2=0.001, f2=0.2)
        with pytest.raises(NotImplementedError):
            KnudsenVia(dual_chem, self.p, self.T, self.AR)

    def test_t0_positive(self):
        """Test that t0() returns a positive characteristic time."""
        model = KnudsenVia(self.chem, self.p, self.T, self.AR)
        t0 = model.t0()
        assert t0 > 0

    def test_saturation_curve(self):
        """Test saturation_curve returns proper data."""
        model = KnudsenVia(self.chem, self.p, self.T, self.AR)
        t, z, cov = model.saturation_curve()
        # Should have multiple time points
        assert len(t) > 1
        assert len(cov) > 1
        # z should have spatial coordinates
        assert len(z) > 1
        # Time should be monotonically increasing
        for i in range(len(t)-1):
            assert t[i] < t[i+1]
        # Coverage arrays should have spatial profiles (each cov element is an array)
        assert hasattr(cov[-1], '__len__') and len(cov[-1]) > 1

    def test_run_default(self):
        """Test run() with default parameters (saturation curve)."""
        model = KnudsenVia(self.chem, self.p, self.T, self.AR)
        t, z, cov = model.run()
        # Should return same as saturation_curve
        assert len(t) > 1
        assert len(z) > 1
        assert len(cov) > 1

    def test_run_with_tdose(self):
        """Test run() with specified dose time."""
        model = KnudsenVia(self.chem, self.p, self.T, self.AR)
        tdose = model.t0() * 2  # Run for 2x characteristic time
        t, z, cov = model.run(tdose=tdose)
        assert len(t) > 1
        assert len(z) > 1
        assert len(cov) > 1

    def test_run_with_tdose_and_dt(self):
        """Test run() with specified dose time and time step."""
        model = KnudsenVia(self.chem, self.p, self.T, self.AR)
        tdose = model.t0() * 2
        dtrun = 0.01
        t, z, cov = model.run(tdose=tdose, dtrun=dtrun)
        assert len(t) > 1
        assert len(z) > 1
        assert len(cov) > 1

    def test_high_aspect_ratio(self):
        """Test model with high aspect ratio."""
        model = KnudsenVia(self.chem, self.p, self.T, AR=50)
        t, z, cov = model.saturation_curve()
        # High AR should still produce valid results
        assert len(t) > 1
        assert len(z) > 1
        assert hasattr(cov[-1], '__len__') and len(cov[-1]) > 1
