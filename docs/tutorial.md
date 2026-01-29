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
3. Running simulations to predict coverage evolution and dose requirements

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
ald = ALDchem(1e2, 1e19, 0.001)
```

The first parameter represents the precursor pressure (in Pa), the second is the number of surface sites per surface area (in square meters), and the third parameter is the sticking probability of the self-limiting process.

This shows a key aspect of aldsim: it works in SI units.

### Using custom units

aldsim lists properunits as a dependency. properunits is a python module that provides a simple way of defining variables
in custom units and transform them into SI units. For instance:

```py
from properunits import Pressure, Temperature

pressure = Pressure(0.1, 'Torr')
temperature = Temperature(200, 'K')

ald = ALDChem(pressure.x, 1e19, temperature.x)
```

Note that in order to extract the SI units, we need to use `pressure.x` and `temperature.x`.

### Relating the density of surface sites to experimental observables

aldsim defines a set of utility functions that help calculate the density of surface sites from observables. 

### Surface Coverage

Surface coverage, typically denoted as θ, represents the fraction of available surface sites that have reacted with the precursor. In an ideal self-limiting process:

- θ = 0 corresponds to a completely unreacted surface
- θ = 1 corresponds to saturation, where all available sites have reacted

### Dose and Exposure

The precursor dose is characterized by the exposure, defined as the product of precursor partial pressure and time:

$$E = P \cdot t$$

The relationship between exposure and surface coverage depends on the reaction mechanism and transport conditions. For a simple first-order Langmuir adsorption model:

$$\frac{d\theta}{dt} = k \cdot P \cdot (1 - \theta)$$

where k is the rate constant for surface reaction.

### Saturation Behavior

A key characteristic of ALD is the saturation curve, which shows how coverage approaches unity with increasing exposure. The shape of this curve depends on:

- **Reaction kinetics**: The intrinsic rate of surface reactions
- **Precursor transport**: How efficiently precursors reach the surface
- **Steric effects**: Site blocking by adsorbed species or ligands

Understanding these factors is essential for optimizing dose times and achieving uniform coatings in practical ALD applications.
