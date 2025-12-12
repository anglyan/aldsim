from aldsim.core.particle import FluidizedBed, RotatingDrum, SpatialPlugFlow, SpatialWellMixed

class TestFluidizedBed:

    def test_saturationcurve(self):
        pfm = FluidizedBed(10)
        x, y = pfm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = FluidizedBed(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


class TestRotatingDrum:

    def test_saturationcurve(self):
        wsm = RotatingDrum(10)
        x, y = wsm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = RotatingDrum(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


class TestSpatialPlugFlow:

    def test_saturationcurve(self):
        wsm = SpatialPlugFlow(10)
        x, y = wsm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = SpatialPlugFlow(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


class TestSpatialWellMixed:

    def test_saturationcurve(self):
        wsm = SpatialWellMixed(10)
        x, y = wsm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = SpatialWellMixed(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape

