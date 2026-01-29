#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim.chem import Precursor, SurfaceKinetics, ALDsoft, ALDideal, ALDchem
import pytest


class TestALDchem:

    def test_init_basic(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, 1e19, 0.001)
        assert ald.npaths == 1
        assert ald.single_path == True
        assert ald.has_rec == False
        assert ald.p_stick1 == pytest.approx(0.001)
        assert ald.dm == 1
        assert ald.f1 == 1
        assert ald.f2 == 0
        assert ald.p_stick2 == 0

    def test_init_with_two_channels(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f2=0.2)
        assert ald.f1 == 0.8
        assert ald.f2 == 0.2
        assert ald.p_stick1 == pytest.approx(0.001)
        assert ald.p_stick2 == pytest.approx(0.0001)
        assert ald.p_rec1 == 0
        assert ald.p_rec0 == 0

    def test_init_with_two_channels_two_f(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f1=0.3, f2=0.2)
        assert ald.f1 == 0.3
        assert ald.f2 == 0.2
        assert ald.p_stick1 == pytest.approx(0.001)
        assert ald.p_stick2 == pytest.approx(0.0001)

    def test_init_with_rec(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f1=0.3, f2=0.2, p_rec1=0.001)
        assert ald.has_rec == True
        assert ald.p_rec1 == pytest.approx(0.001)
        assert ald.p_rec0 == 0

    def test_init_with_rec2(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f1=0.3, f2=0.2, p_rec1=0.001, p_rec0=0.01)
        assert ald.has_rec == True
        assert ald.p_rec1 == pytest.approx(0.001)
        assert ald.p_rec0 == pytest.approx(0.01)

    def test_sticking_prob(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001, f1=0.3)
        
        with pytest.raises(ValueError):
            ald.sticking_prob(0.001, 0.1)
        assert ald.sticking_prob(1) == 0
        assert ald.sticking_prob(0) == pytest.approx(0.0003)

    def test_sticking_prob_nof(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001)
        
        assert ald.sticking_prob(0) == pytest.approx(0.001)


    def test_sticking_prob2(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f1=0.3, f2=0.2)
        
        with pytest.raises(ValueError):
            ald.sticking_prob(0.001)
        
        assert ald.sticking_prob(0,0) == pytest.approx(0.00032)


    def test_sticking_prob_av(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001, f1=0.3)
        
        with pytest.raises(ValueError):
            ald.sticking_prob_av(0.001, 0.1)
        assert ald.sticking_prob_av(0) == 0
        assert ald.sticking_prob_av(1) == pytest.approx(0.0003)

    def test_sticking_prob_av_nof(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001)
        
        assert ald.sticking_prob_av(1) == pytest.approx(0.001)


    def test_sticking_prob_av2(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f1=0.3, f2=0.2)
        
        with pytest.raises(ValueError):
            ald.sticking_prob_av(0.001)
        
        assert ald.sticking_prob_av(1,1) == pytest.approx(0.00032)

    def test_recomb_prob(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001, f1=0.3,
            p_rec0=0.01, p_rec1=0.1)

        with pytest.raises(ValueError):
            ald.recomb_prob(0.001, 0.1)
        assert ald.recomb_prob(0) == pytest.approx(0.01)
        assert ald.recomb_prob(1) == pytest.approx(0.037)

    def test_recomb_prob_nof(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_rec0=0.01, p_rec1=0.1)

        assert ald.recomb_prob(0) == pytest.approx(0.01)
        assert ald.recomb_prob(1) == pytest.approx(0.1)

    def test_recomb_prob2(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f1=0.3, f2=0.2, p_rec0=0.01, p_rec1=0.1)

        with pytest.raises(ValueError):
            ald.recomb_prob(0.001)

        assert ald.recomb_prob(0, 0) == pytest.approx(0.01)
        assert ald.recomb_prob(1, 1) == pytest.approx(0.055)

    def test_react_prob(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001, f1=0.3,
            p_rec0=0.01, p_rec1=0.1)

        with pytest.raises(ValueError):
            ald.react_prob(0.001, 0.1)
        assert ald.react_prob(0) == pytest.approx(0.0103)
        assert ald.react_prob(1) == pytest.approx(0.037)

    def test_react_prob_nof(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_rec0=0.01, p_rec1=0.1)

        assert ald.react_prob(0) == pytest.approx(0.011)
        assert ald.react_prob(1) == pytest.approx(0.1)

    def test_react_prob2(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f1=0.3, f2=0.2, p_rec0=0.01, p_rec1=0.1)

        with pytest.raises(ValueError):
            ald.react_prob(0.001)

        assert ald.react_prob(0, 0) == pytest.approx(0.01032)
        assert ald.react_prob(1, 1) == pytest.approx(0.055)


    def test_t0(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001, f1=0.3)
        assert isinstance(ald.t0(200, 1e2), float)

    def test_t02(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001,
            p_stick2=0.0001, f1=0.3, f2=0.2)

        assert isinstance(ald.t0(200, 1e2), tuple)

    def test_saturation(self):
        p = Precursor(mass=100)
        ald = ALDchem(p, nsites=1e19, p_stick1=0.001, f1=0.3)
        assert isinstance(ald.saturation_curve(400, 0.1), tuple)

