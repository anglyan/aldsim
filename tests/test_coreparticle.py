#Copyright © 2024-Present, UChicago Argonne, LLC

from aldsim.core.particle import FluidizedBedND, WellMixedParticleND, SpatialPlugFlow, SpatialWellMixedND

class TestFluidizedBedND:

    def test_saturationcurve(self):
        pfm = FluidizedBedND(10)
        x, y = pfm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = FluidizedBedND(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


class TestWellMixedParticleND:

    def test_saturationcurve(self):
        wsm = WellMixedParticleND(10)
        x, y = wsm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = WellMixedParticleND(10)
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
        wsm = SpatialWellMixedND(10)
        x, y = wsm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = SpatialWellMixedND(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape

