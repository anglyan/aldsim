Core models
===========

aldsim implements a number of core models for atomic layer deposition.

Particle coating
----------------

aldsim implements four particle coating models for atomic layer deposition (ALD) processes, each representing different reactor configurations and transport regimes:

**FluidizedBed**
    Batch fluidized bed reactor model combining well-mixed particle approximation with plug flow precursor transport. Suitable for modeling particle coating in fluidized bed reactors where particles are well mixed but precursor flow exhibits plug flow behavior.

**RotatingDrum**
    Batch particle coating model under well-stirred reactor approximation for both particles and precursor transport. Appropriate for rotating drum reactors and other systems where complete mixing of both phases occurs.

**SpatialPlugFlow**
    Continuous particle coating model with stratified particle mixing and plug flow precursor transport. Both particles and precursor move in the same direction, making this suitable for spatial ALD systems with co-current flow.

**SpatialWellMixed**
    Continuous particle coating model with stratified particle mixing and well-stirred precursor transport. This model inherits from RotatingDrum and is formally equivalent to batch coating when residence time replaces dose time.

All models assume first-order irreversible Langmuir kinetics with reaction rates characterized by a Damkohler number (Da). Each class provides methods to calculate surface coverage, generate saturation curves, and run full simulations over normalized time or residence time.

These models are described in the paper: *Modeling scale up of particle coating by atomic layer deposition* (https://doi.org/10.1116/6.0004006)
by Angel Yanguas-Gil and Jeffrey W. Elam, Journal of Vacuum Science and Technology A 43, 012404 (2025). Please cite this work if
you use any of these models in your research.

Example
-------

This is an example of how to use the `FluidizedBed` class to
compute saturation curve for a given value of the Damkohler number:

.. code-block:: python
    import matplotlib.pyplot as plt
    from aldsim.core.particle import FluidizedBed

    # Create a FluidizedBed instance with Damkohler number = 1
    fb = FluidizedBed(Da=1.0)

    # Calculate the saturation curve
    t, coverage = fb.saturation_curve(tmax=5, dt=0.01)

    # Plot the saturation curve
    plt.figure(figsize=(8, 6))
    plt.plot(t, coverage, 'b-', linewidth=2, label='Da = 1.0')
    plt.xlabel('Normalized Dose Time', fontsize=12)
    plt.ylabel('Surface Coverage', fontsize=12)
    plt.title('Saturation Curve for Fluidized Bed ALD', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
