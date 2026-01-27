#Copyright © 2024, UChicago Argonne, LLC

from .chem import ALDideal
from .models import ZeroD, RotatingDrum, ParticlePlugFlow

_ideal_models = {
    'zeroD' : ZeroD,
    'rotatingdrum' : RotatingDrum,
    'fluidizedbed' : ParticlePlugFlow
}

def aldmodel(process, model_name, **kwargs):
    if isinstance(process, ALDideal):
        return _ideal_models[model_name](process, **kwargs)

