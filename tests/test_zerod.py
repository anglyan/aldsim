#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim import Precursor, aldmodel
from aldsim.chem import ALDideal
from aldsim.models import ZeroD
from aldsim.constants import kb
import numpy as np
import pytest


class TestZeroD:

    def test_init_basic(self):
        """Test basic initialization of ZeroD model"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)
        model = ZeroD(chem, T=500, p=0.1*1e5/760)

        assert model.chem == chem
        assert model.T == 500
        assert model.p == pytest.approx(0.1*1e5/760)

    def test_init_via_aldmodel(self):
        """Test initialization of ZeroD via aldmodel factory function"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)
        model = aldmodel(chem, 'zeroD', p=0.1*1e5/760, T=500)

        assert isinstance(model, ZeroD)
        assert model.T == 500
        assert model.p == pytest.approx(0.1*1e5/760)

    def test_properties_inherited_from_base(self):
        """Test that ZeroD inherits properties from IdealDoseModel"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)
        model = ZeroD(chem, T=500, p=0.1*1e5/760)

        # Test that vth property works
        assert hasattr(model, 'vth')
        assert model.vth > 0

    def test_saturation_curve_returns_arrays(self):
        """Test that saturation_curve returns proper numpy arrays"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)
        model = ZeroD(chem, T=500, p=0.1*1e5/760)

        t_arr, cov_arr = model.saturation_curve()

        # Check that both are numpy arrays
        assert isinstance(t_arr, np.ndarray)
        assert isinstance(cov_arr, np.ndarray)

        # Check that they have the same length
        assert len(t_arr) == len(cov_arr)

    def test_saturation_curve_time_properties(self):
        """Test time array properties from saturation_curve"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)
        model = ZeroD(chem, T=500, p=0.1*1e5/760)

        t_arr, cov_arr = model.saturation_curve()

        # Time should start at 0
        assert t_arr[0] == 0

        # Time should be monotonically increasing
        assert all(t_arr[i] < t_arr[i+1] for i in range(len(t_arr)-1))

        # All times should be non-negative
        assert all(t >= 0 for t in t_arr)

    def test_saturation_curve_coverage_properties(self):
        """Test coverage array properties from saturation_curve"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)
        model = ZeroD(chem, T=500, p=0.1*1e5/760)

        t_arr, cov_arr = model.saturation_curve()

        # Coverage should start at 0
        assert cov_arr[0] == pytest.approx(0)

        # Coverage should be monotonically increasing
        assert all(cov_arr[i] <= cov_arr[i+1] for i in range(len(cov_arr)-1))

        # Coverage should be bounded between 0 and 1
        assert all(0 <= c <= 1 for c in cov_arr)

        # Coverage should approach 1 asymptotically (should be > 0.9 at end)
        assert cov_arr[-1] > 0.9

    def test_saturation_curve_exponential_behavior(self):
        """Test that saturation curve follows expected exponential form"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)
        model = ZeroD(chem, T=500, p=0.1*1e5/760)

        t_arr, cov_arr = model.saturation_curve()

        # Calculate expected behavior: cov = 1 - exp(-t/t0)
        # t0 = 1/nu where nu is calculated in the method
        nu = 0.25 * model.chem.site_area * model.vth * model.p / (kb * model.T) * model.chem.p_stick0
        t0 = 1/nu
        expected_cov = 1 - np.exp(-t_arr/t0)

        # Check that calculated coverage matches expected exponential
        assert np.allclose(cov_arr, expected_cov, rtol=1e-10)

    def test_saturation_curve_different_temperatures(self):
        """Test saturation curve at different temperatures"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)

        # Lower temperature
        model_low = ZeroD(chem, T=400, p=0.1*1e5/760)
        t_low, cov_low = model_low.saturation_curve()

        # Higher temperature
        model_high = ZeroD(chem, T=600, p=0.1*1e5/760)
        t_high, cov_high = model_high.saturation_curve()

        # Both should have proper coverage behavior
        assert cov_low[-1] > 0.9
        assert cov_high[-1] > 0.9

    def test_saturation_curve_different_pressures(self):
        """Test saturation curve at different pressures"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)

        # Lower pressure
        model_low = ZeroD(chem, T=500, p=0.05*1e5/760)
        t_low, cov_low = model_low.saturation_curve()

        # Higher pressure
        model_high = ZeroD(chem, T=500, p=0.2*1e5/760)
        t_high, cov_high = model_high.saturation_curve()

        # Both should saturate
        assert cov_low[-1] > 0.9
        assert cov_high[-1] > 0.9

        # Higher pressure should have shorter characteristic time
        # (saturation happens faster)
        nu_low = 0.25 * model_low.chem.site_area * model_low.vth * model_low.p / (kb * model_low.T) * model_low.chem.p_stick0
        nu_high = 0.25 * model_high.chem.site_area * model_high.vth * model_high.p / (kb * model_high.T) * model_high.chem.p_stick0
        assert nu_high > nu_low  # Higher pressure -> faster saturation rate

    def test_saturation_curve_different_sticking_coefficients(self):
        """Test saturation curve with different sticking coefficients"""
        prec = Precursor(mass=150.0)

        # Lower sticking coefficient
        chem_low = ALDideal(prec, 1e19, 1e-4)
        model_low = ZeroD(chem_low, T=500, p=0.1*1e5/760)
        t_low, cov_low = model_low.saturation_curve()

        # Higher sticking coefficient
        chem_high = ALDideal(prec, 1e19, 1e-2)
        model_high = ZeroD(chem_high, T=500, p=0.1*1e5/760)
        t_high, cov_high = model_high.saturation_curve()

        # Both should saturate
        assert cov_low[-1] > 0.9
        assert cov_high[-1] > 0.9

        # Higher sticking coefficient should have faster saturation
        nu_low = 0.25 * model_low.chem.site_area * model_low.vth * model_low.p / (kb * model_low.T) * model_low.chem.p_stick0
        nu_high = 0.25 * model_high.chem.site_area * model_high.vth * model_high.p / (kb * model_high.T) * model_high.chem.p_stick0
        assert nu_high > nu_low

    def test_temperature_change_updates_vth(self):
        """Test that changing temperature updates thermal velocity"""
        prec = Precursor(mass=150.0)
        chem = ALDideal(prec, 1e19, 1e-3)
        model = ZeroD(chem, T=500, p=0.1*1e5/760)

        vth_initial = model.vth

        # Change temperature
        model.T = 600

        # vth should change
        assert model.vth != vth_initial
        assert model.vth > vth_initial  # Higher T -> higher vth

    def test_basic_example_workflow(self):
        """Test the workflow from basic_example.py"""
        prec = Precursor(mass=150.0)
        nsites = 1e19
        beta0 = 1e-3
        chem = ALDideal(prec, nsites, beta0, dm=1.0)
        model = aldmodel(chem, 'zeroD', p=0.1*1e5/760, T=500)

        t, theta = model.saturation_curve()

        # Verify the output is reasonable
        assert len(t) > 0
        assert len(theta) > 0
        assert len(t) == len(theta)
        assert theta[-1] > 0.9
        assert theta[0] == pytest.approx(0)
