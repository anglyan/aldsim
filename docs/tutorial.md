# Tutorial

## Introduction

This tutorial provides a practical guide to using aldsim for modeling atomic layer deposition processes.

Before starting, ensure aldsim is installed:

```python
pip install aldsim
```

The basic workflow in aldsim involves:

1. Defining the surface chemistry and kinetic parameters
2. Selecting an appropriate reactor or transport model
3. Running simulations to model the evolution of surface coverage with time

## ALD Surface kinetics

The kinetics of ALD processes are governed by the interaction between precursor molecules and reactive surface sites. In aldsim, these interactions are modeled through the `chem` module, which provides abstractions for common self-limiting reaction mechanisms.

The first step in defining a specific ALD process is
to create an ALD precursor as follows:

```py
import aldsim

prec = aldsim.Precursor(mass=100)
```

The parameter `mass` provides the precursor mass in atomic
mass units. 

We can now define a self-limited process through the class `ALDchem`:

```py
ald = ALDchem(prec, 1e19, 0.001)
```

The first parameter represents the precursor, the second is the number of surface sites per surface area (in square meters), and the third parameter is the sticking probability of the self-limiting process.

### Relating the density of surface sites to experimental observables

aldsim defines a set of utility functions that help calculate the density of surface sites from observables. 
The first function is `sitedensity_fromqcm`, which calculates the density of surface sites from the mass (in ng/cm2)
obtained by quartz crystal microbalance. It requires two additional variables: the molar mass of the solid (in g)
and the number of precursor molecules required to complete the molecular formula. The idea here is that,
given a mass per surface area and  the mass per moles, you can get the moles per surface area and therefore
the number of sites per surface area. For instance, for TMA and alumina:

```py
from aldsim.utils import sitedensity_fromqcm

mm_Al2O3 = 102
s1 = sitedensity_fromqcm(35, 102, N=2)
```

The value of `s1` is approximately `4e18` sites/m2.

## A simple reactor model

Once we have defined a surface kinetics we can load an ALD model to compute the evolution of surface coverage as a function
of time. Here is an example for the case of a 0D model that models the evolution of surface coverage as a function of time
for a constant precursor partial pressure:

```py
from aldsim.models import ZeroD
model = ZeroD(ald, T=500, p=1e2)
t, cov = model.saturation_curve()
```
Note that temperature and pressure both use SI units, so the temperature is in K and the pressure is in Pa.
The method `saturation_curve` returns the saturation curve for the self-limited process for these two conditions as
a tuple of two `numpy.array` objects.

### Using custom units

aldsim lists properunits as a dependency. properunits is a python module that provides a simple way of defining variables
in custom units and transform them into SI units. For instance:

```py
from properunits import Pressure, Temperature

pressure = Pressure(0.1, 'Torr')
temperature = Temperature(200, 'C')
model = ZeroD(ald, T=temperature.x, p=pressure.x)
```

In order to extract the SI units, we need to use `pressure.x` and `temperature.x`. This provides a simple way of
defining these variables in other unitt, such as degree celsius or Torr.

## Example: particle coating by ALD

Here we include a simple example of how to model ALD in a rotating drum reactor under a well mixed approximation.

```py
from aldsim import Precursor, ALDchem
from aldsim.models import RotatingDrum

import matplotlib.pyplot as pt

prec = Precursor(mass=150.0)
nsites = 1e19
beta0 = 1e-3
chem = ALDchem(prec, nsites, beta0, dm=1.0)
model = RotatingDrum(chem, p=0.1*1e5/760, p0=1e2, T=500,
                 S=1e1, flow=60)

t, theta = model.saturation_curve()

pt.plot(t, theta)
pt.show()
```

In this example, `RotatingDrum` takes an ALDchem object, the precursor pressure, the base pressure, temperature,
the total surface area (in square meters) and the flow to the reactor in sccm. 
