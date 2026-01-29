#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim.core.particle import SpatialPlugFlow
from aldsim.core.particle.spatialplugflow import saturation_curve
import numpy as np
import pytest


class TestSpatialPlugFlow:

    def test_init_basic(self):
        """Test basic initialization of SpatialPlugFlow"""
        model = SpatialPlugFlow(Da=1.0)
        assert model.Da == 1.0

    def test_init_different_damkohler(self):
        """Test initialization with different Damkohler numbers"""
        model1 = SpatialPlugFlow(Da=0.5)
        model2 = SpatialPlugFlow(Da=10.0)
        assert model1.Da == 0.5
        assert model2.Da == 10.0

    # calc_coverage tests
    def test_calc_coverage_at_unity_time(self):
        """Test calc_coverage at t=1 (residence time equals reactor time)"""
        model = SpatialPlugFlow(Da=2.0)
        coverage = model.calc_coverage(t=1)
        # At t=1, coverage = Da/(1+Da)
        expected = 2.0 / (1 + 2.0)
        assert coverage == pytest.approx(expected)

    def test_calc_coverage_at_unity_time_various_da(self):
        """Test calc_coverage at t=1 for various Damkohler numbers"""
        da_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
        for da in da_values:
            model = SpatialPlugFlow(Da=da)
            coverage = model.calc_coverage(t=1)
            expected = da / (1 + da)
            assert coverage == pytest.approx(expected)

    def test_calc_coverage_at_non_unity_time(self):
        """Test calc_coverage at t != 1"""
        model = SpatialPlugFlow(Da=2.0)
        coverage = model.calc_coverage(t=0.5)
        # At t != 1, coverage = 1 - (1-t)/(1-t*exp(-Da*(1-t)))
        x = np.exp(-2.0 * (1 - 0.5))
        expected = 1 - (1 - 0.5) / (1 - 0.5 * x)
        assert coverage == pytest.approx(expected)

    def test_calc_coverage_time_zero(self):
        """Test calc_coverage at t=0 (no residence time)"""
        model = SpatialPlugFlow(Da=2.0)
        coverage = model.calc_coverage(t=0)
        # At t=0, should give zero coverage
        x = np.exp(-2.0 * 1)
        expected = 1 - 1 / (1 - 0)
        assert coverage == pytest.approx(0.0)

    def test_calc_coverage_increases_with_time(self):
        """Test that coverage increases monotonically with time"""
        model = SpatialPlugFlow(Da=2.0)
        t_values = np.linspace(0.1, 3.0, 20)
        coverages = [model.calc_coverage(t=t) for t in t_values]
        # Coverage should be monotonically increasing
        assert all(coverages[i] <= coverages[i+1] for i in range(len(coverages)-1))

    def test_calc_coverage_bounded(self):
        """Test that coverage is bounded between 0 and 1"""
        model = SpatialPlugFlow(Da=2.0)
        t_values = np.linspace(0.01, 10.0, 100)
        coverages = [model.calc_coverage(t=t) for t in t_values]
        # All coverages should be in [0, 1]
        assert all(0 <= c <= 1 for c in coverages)

    def test_calc_coverage_with_da_override(self):
        """Test calc_coverage with Da parameter override"""
        model = SpatialPlugFlow(Da=1.0)
        coverage = model.calc_coverage(t=1, Da=2.0)
        # Should use the overridden Da=2.0
        expected = 2.0 / (1 + 2.0)
        assert coverage == pytest.approx(expected)
        # Model's Da should be updated
        assert model.Da == 2.0

    def test_calc_coverage_approaches_saturation(self):
        """Test that coverage approaches saturation at large times"""
        model = SpatialPlugFlow(Da=2.0)
        cov_large_t = model.calc_coverage(t=10.0)
        # Should be close to 1
        assert cov_large_t > 0.95

    # calc_precursor tests
    def test_calc_precursor_at_unity_time(self):
        """Test calc_precursor at t=1"""
        model = SpatialPlugFlow(Da=2.0)
        prec_unused = model.calc_precursor(t=1)
        # At t=1, precursor unused fraction = 1 - Da/(1+Da)
        expected = 1 - 2.0 / (1 + 2.0)
        assert prec_unused == pytest.approx(expected)

    def test_calc_precursor_at_non_unity_time(self):
        """Test calc_precursor at t != 1"""
        model = SpatialPlugFlow(Da=2.0)
        prec_unused = model.calc_precursor(t=0.5)
        # At t != 1, precursor unused = (1-t)*exp(-Da*(1-t))/(1-t*exp(-Da*(1-t)))
        x = np.exp(-2.0 * (1 - 0.5))
        expected = (1 - 0.5) * x / (1 - 0.5 * x)
        assert prec_unused == pytest.approx(expected)

    def test_calc_precursor_bounded(self):
        """Test that precursor unused fraction is bounded between 0 and 1"""
        model = SpatialPlugFlow(Da=2.0)
        t_values = np.linspace(0.01, 5.0, 100)
        prec_values = [model.calc_precursor(t=t) for t in t_values]
        # All precursor values should be in [0, 1]
        assert all(0 <= p <= 1 for p in prec_values)

    def test_calc_precursor_with_da_override(self):
        """Test calc_precursor with Da parameter override"""
        model = SpatialPlugFlow(Da=1.0)
        prec_unused = model.calc_precursor(t=1, Da=2.0)
        # Should use the overridden Da=2.0
        expected = 1 - 2.0 / (1 + 2.0)
        assert prec_unused == pytest.approx(expected)
        # Model's Da should be updated
        assert model.Da == 2.0

    def test_precursor_coverage_conservation(self):
        """Test that coverage and precursor usage are related properly"""
        model = SpatialPlugFlow(Da=2.0)
        t = 1.0
        coverage = model.calc_coverage(t=t)
        prec_unused = model.calc_precursor(t=t)
        # At t=1, coverage + prec_unused should equal 1
        assert coverage + prec_unused == pytest.approx(1.0)

    # saturation_curve tests
    def test_saturation_curve_returns_arrays(self):
        """Test that saturation_curve returns numpy arrays"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage = model.saturation_curve()
        assert isinstance(t, np.ndarray)
        assert isinstance(coverage, np.ndarray)

    def test_saturation_curve_same_length(self):
        """Test that time and coverage arrays have same length"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage = model.saturation_curve()
        assert len(t) == len(coverage)

    def test_saturation_curve_default_parameters(self):
        """Test saturation_curve with default parameters"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage = model.saturation_curve()
        # Default: tmax=5, dt=0.01
        assert t[0] == 0.0
        assert t[-1] == pytest.approx(5.0 - 0.01)
        assert len(t) == int(5.0 / 0.01)

    def test_saturation_curve_custom_parameters(self):
        """Test saturation_curve with custom parameters"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage = model.saturation_curve(tmax=3.0, dt=0.05)
        assert t[0] == 0.0
        assert t[-1] == pytest.approx(3.0 - 0.05)
        assert len(t) == int(3.0 / 0.05)

    def test_saturation_curve_time_properties(self):
        """Test time array properties from saturation_curve"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage = model.saturation_curve()
        # Time should start at 0
        assert t[0] == 0.0
        # Time should be monotonically increasing
        assert all(t[i] < t[i+1] for i in range(len(t)-1))
        # All times should be non-negative
        assert all(ti >= 0 for ti in t)

    def test_saturation_curve_coverage_properties(self):
        """Test coverage array properties from saturation_curve"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage = model.saturation_curve()
        # Coverage should be monotonically increasing
        assert all(coverage[i] <= coverage[i+1] for i in range(len(coverage)-1))
        # Coverage should be bounded between 0 and 1
        assert all(0 <= c <= 1 for c in coverage)
        # Coverage should start near 0
        assert coverage[0] < 0.1

    def test_saturation_curve_different_damkohler(self):
        """Test saturation curves for different Damkohler numbers"""
        model_low = SpatialPlugFlow(Da=0.5)
        model_high = SpatialPlugFlow(Da=5.0)

        t_low, cov_low = model_low.saturation_curve()
        t_high, cov_high = model_high.saturation_curve()

        # At same residence time, higher Da should give higher coverage
        for i in range(len(t_low)):
            if t_low[i] > 0:  # Skip t=0 where both are 0
                assert cov_high[i] >= cov_low[i]

    def test_saturation_curve_module_function(self):
        """Test the module-level saturation_curve function"""
        t, coverage = saturation_curve(Da=2.0)
        assert isinstance(t, np.ndarray)
        assert isinstance(coverage, np.ndarray)
        assert len(t) == len(coverage)

    def test_saturation_curve_module_function_custom_params(self):
        """Test module-level saturation_curve with custom parameters"""
        t, coverage = saturation_curve(Da=2.0, tmax=3.0, dt=0.1)
        assert t[0] == 0.0
        assert t[-1] == pytest.approx(3.0 - 0.1)

    # run tests
    def test_run_returns_three_arrays(self):
        """Test that run returns three arrays"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage, precursor = model.run()
        assert isinstance(t, np.ndarray)
        assert isinstance(coverage, np.ndarray)
        assert isinstance(precursor, np.ndarray)

    def test_run_same_length_arrays(self):
        """Test that all arrays from run have same length"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage, precursor = model.run()
        assert len(t) == len(coverage)
        assert len(t) == len(precursor)

    def test_run_default_parameters(self):
        """Test run with default parameters"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage, precursor = model.run()
        # Default: tmax=5, dt=0.01
        assert t[0] == 0.0
        assert t[-1] == pytest.approx(5.0 - 0.01)
        assert len(t) == int(5.0 / 0.01)

    def test_run_custom_parameters(self):
        """Test run with custom parameters"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage, precursor = model.run(tmax=3.0, dt=0.05)
        assert t[0] == 0.0
        assert t[-1] == pytest.approx(3.0 - 0.05)
        assert len(t) == int(3.0 / 0.05)

    def test_run_time_properties(self):
        """Test time array properties from run"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage, precursor = model.run()
        # Time should start at 0
        assert t[0] == 0.0
        # Time should be monotonically increasing
        assert all(t[i] < t[i+1] for i in range(len(t)-1))
        # All times should be non-negative
        assert all(ti >= 0 for ti in t)

    def test_run_coverage_properties(self):
        """Test coverage array properties from run"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage, precursor = model.run()
        # Coverage should be monotonically increasing
        assert all(coverage[i] <= coverage[i+1] for i in range(len(coverage)-1))
        # Coverage should be bounded between 0 and 1
        assert all(0 <= c <= 1 for c in coverage)


    def test_run_precursor_coverage_relationship(self):
        """Test relationship between coverage and precursor utilization from run"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage, precursor = model.run()

        # For each point where t is close to 1, coverage + precursor should be ~1
        tolerance = 0.01
        for i, ti in enumerate(t):
            if abs(ti - 1.0) < tolerance:
                assert coverage[i] + precursor[i] == pytest.approx(1.0, abs=0.01)
                break

    def test_run_different_damkohler(self):
        """Test run for different Damkohler numbers"""
        model_low = SpatialPlugFlow(Da=0.5)
        model_high = SpatialPlugFlow(Da=5.0)

        t_low, cov_low, prec_low = model_low.run()
        t_high, cov_high, prec_high = model_high.run()

        # At same residence time, higher Da should give:
        # - Higher coverage
        # - Lower unused precursor (more utilized)
        for i in range(len(t_low)):
            if t_low[i] > 0:  # Skip t=0
                assert cov_high[i] >= cov_low[i]
                assert prec_high[i] <= prec_low[i]

    def test_run_matches_individual_calculations(self):
        """Test that run results match individual calc_coverage and calc_precursor calls"""
        model = SpatialPlugFlow(Da=2.0)
        t, coverage, precursor = model.run(tmax=2.0, dt=0.1)

        # Check a few random points
        test_indices = [0, 5, 10, 15, 19]
        for i in test_indices:
            ti = t[i]
            expected_cov = model.calc_coverage(t=ti)
            expected_prec = model.calc_precursor(t=ti)
            assert coverage[i] == pytest.approx(expected_cov)
            assert precursor[i] == pytest.approx(expected_prec)

    def test_run_consistency_with_saturation_curve(self):
        """Test that run and saturation_curve give consistent coverage results"""
        model = SpatialPlugFlow(Da=2.0)

        t_run, cov_run, prec_run = model.run(tmax=3.0, dt=0.02)
        t_sat, cov_sat = model.saturation_curve(tmax=3.0, dt=0.02)

        # Time arrays should be identical
        assert np.allclose(t_run, t_sat)
        # Coverage arrays should be identical
        assert np.allclose(cov_run, cov_sat)
