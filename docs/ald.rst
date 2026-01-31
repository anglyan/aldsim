About Atomic Layer Deposition
=============================

Atomic layer deposition (ALD) is a thin film growth technique based on sequential, self-limiting surface reactions. ALD achieves precise thickness control at the atomic scale through cyclic exposure of precursor gases.

Anatomy of a typical ALD Process
--------------------------------

A typical ALD cycle consists of four steps:

1. **Precursor dose**: The first precursor is introduced and adsorbs onto the surface through chemisorption
2. **Purge**: Excess precursor and reaction byproducts are removed
3. **Co-reactant dose**: A second precursor reacts with the adsorbed species
4. **Purge**: The reactor is purged again to complete the cycle

The self-limiting nature of each half-reaction ensures that film growth is controlled by the number of cycles rather than process parameters like temperature or pressure, enabling conformal coating on complex geometries.

aldsim focuses primarily on modeling ALD processes during either a precursor or a coreactant dose.

ALD surface kinetics in aldsim
------------------------------

aldsim implements self-limited surface kinetics through its ``ALDchem`` class. ``ALDchem`` implements surface kinetics with up
to four reaction pathways: two self-limited reaction pathways used to model soft-saturating processes, and two recombination
pathways depending on whether a surface site is available or occupied by a dissociatively adsorbed precursor. These
recombination pathways are used to model plasma and ozone-based ALD processes.
