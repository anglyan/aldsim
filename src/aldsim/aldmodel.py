#Copyright © 2024-Present, UChicago Argonne, LLC

from .chem import ALDideal
from .models import ZeroD, RotatingDrum, FluidizedBed

_ideal_models = {
    'zeroD' : ZeroD,
    'rotatingdrum' : RotatingDrum,
    'fluidizedbed' : FluidizedBed
}

def aldmodel(process, model_name, **kwargs):
    if isinstance(process, ALDideal):
        return _ideal_models[model_name](process, **kwargs)

