from aldsim.chem import Precursor, SurfaceKinetics, ALDsoft, ALDideal
import pytest


def test_precursor():
    c = Precursor(mass=100)
    assert c.name == 'None'
    assert c.mass == 100
    c = Precursor(name='TMA', mass=101.5)
    assert c.name == 'TMA'
        

class TestSurfaceKinetics:

    def test_surfacekinetics(self):
        p = Precursor(mass=100)
        k = SurfaceKinetics(p, 1e19, 1)
        assert k.site_area == pytest.approx(1e-19)

    def test_surfacekinupdate(self):
        p = Precursor(mass=100)
        k = SurfaceKinetics(p, 1e19, 1)
        k.site_area = 1e-18
        assert k.nsites == pytest.approx(1e18)


class TestALDideal:

    def test_init_basic(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001)
        assert ald.p_stick0 == 0.001
        assert ald.dm == 1
        assert ald.f == 1

    def test_init_with_custom_f(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001, f=0.8)
        assert ald.f == 0.8
        assert ald.p_stick0 == 0.001

    def test_init_with_custom_dm(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001, dm=2)
        assert ald.dm == 2

    def test_site_area(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001)
        assert ald.site_area == pytest.approx(1e-19)

    def test_sticking_prob_zero_coverage(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001)
        # At zero coverage, sprob = f*sprob0
        expected = 1.0 * 0.001
        assert ald.sticking_prob(0) == pytest.approx(expected)

    def test_sticking_prob_partial_coverage(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001)
        # At partial coverage cov=0.5
        expected = 1.0 * 0.001 * (1-0.5)
        assert ald.sticking_prob(0.5) == pytest.approx(expected)

    def test_sticking_prob_full_coverage(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001)
        # At full coverage, sprob should be zero
        assert ald.sticking_prob(1) == pytest.approx(0)

    def test_sticking_prob_with_fractional_f(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001, f=0.8)
        # At zero coverage with f=0.8
        expected = 0.8 * 0.001
        assert ald.sticking_prob(0) == pytest.approx(expected)

    def test_sticking_prob_av(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001)
        # Test sprob_av with average value
        av = 0.6
        expected = 1.0 * 0.001 * av
        assert ald.sticking_prob_av(av) == pytest.approx(expected)

    def test_t0(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001)
        T = 300  # K
        pressure = 100  # Pa
        t0 = ald.t0(T, pressure)
        # Check that t0 is inversely proportional to sprob0
        expected_t0 = 1.0/(ald.site_area * ald.Jwall(T, pressure) * ald.p_stick0)
        assert t0 == pytest.approx(expected_t0)
        # t0 should be positive
        assert t0 > 0

    def test_saturation_curve(self):
        p = Precursor(mass=100)
        ald = ALDideal(p, 1e19, 0.001)
        T = 300  # K
        pressure = 100  # Pa
        t, cov = ald.saturation_curve(T, pressure)
        # Check that time array is monotonically increasing
        assert all(t[i] < t[i+1] for i in range(len(t)-1))
        # Check that coverage starts near zero
        assert cov[0] < 0.1
        # Check that coverage is bounded between 0 and 1
        assert all(0 <= c <= 1 for c in cov)
        # Check that coverage increases monotonically
        assert all(cov[i] <= cov[i+1] for i in range(len(cov)-1))
        # For ideal kinetics, should approach 1 asymptotically
        assert cov[-1] > 0.9


class TestSoftSaturating:

    def test_nsites_s0(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        assert k.site_area == pytest.approx(1e-19)
        k.site_area = 1e-18
        assert k.nsites == pytest.approx(1e18)

    def test_init_with_explicit_f2(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        assert k.f1 == 0.8
        assert k.f2 == 0.2
        assert k.p_stick1 == 1e-2
        assert k.p_stick2 == 1e-3

    def test_init_with_implicit_f2(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8)
        assert k.f1 == 0.8
        assert k.f2 == pytest.approx(0.2)
        assert k.f == pytest.approx(1.0)

    def test_sticking_prob_zero_coverage(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        # At zero coverage, sticking_prob = f1*p_stick1 + f2*p_stick2
        expected = 0.8 * 1e-2 + 0.2 * 1e-3
        assert k.sticking_prob(0, 0) == pytest.approx(expected)

    def test_sticking_prob_partial_coverage(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        # At partial coverage cov1=0.5, cov2=0.3
        expected = 0.8 * 1e-2 * (1-0.5) + 0.2 * 1e-3 * (1-0.3)
        assert k.sticking_prob(0.5, 0.3) == pytest.approx(expected)

    def test_sticking_prob_full_coverage(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        # At full coverage, sticking_prob should be zero
        assert k.sticking_prob(1, 1) == pytest.approx(0)

    def test_sticking_prob_av(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        # Test sticking_prob_av with average values
        av1, av2 = 0.6, 0.4
        expected = 0.8 * 1e-2 * av1 + 0.2 * 1e-3 * av2
        assert k.sticking_prob_av(av1, av2) == pytest.approx(expected)

    def test_t0(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        T = 300  # K
        pressure = 100  # Pa
        t1, t2 = k.t0(T, pressure)
        # t1 should be smaller than t2 (since p_stick1 > p_stick2)
        assert t1 < t2
        # Check that t1 is inversely proportional to p_stick1
        expected_t1 = 1.0/(k.site_area * k.Jwall(T, pressure) * k.p_stick1)
        expected_t2 = 1.0/(k.site_area * k.Jwall(T, pressure) * k.p_stick2)
        assert t1 == pytest.approx(expected_t1)
        assert t2 == pytest.approx(expected_t2)

    def test_saturation_curve(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        T = 300  # K
        pressure = 100  # Pa
        t, cov = k.saturation_curve(T, pressure)
        # Check that time array is monotonically increasing
        assert all(t[i] < t[i+1] for i in range(len(t)-1))
        # Check that coverage starts near zero and increases
        assert cov[0] < 0.1
        # Check that coverage is bounded between 0 and 1
        assert all(0 <= c <= 1 for c in cov)
        # Check that coverage increases monotonically
        assert all(cov[i] <= cov[i+1] for i in range(len(cov)-1))
