.. aldsim documentation master file, created by
   sphinx-quickstart on Sat May  4 15:42:36 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

aldsim's documentation
======================

aldsim provides a suite of simple models for atomic layer deposition
processes.

Atomic layer deposition (ALD) is a thin film growth technique that
relies on self-limited surface kinetics. It plays a key role
in areas such as microelectronics, and it is applied for energy,
energy storage, catalysis, and decarbonization applications.

aldsim implements a series of models to help explore ALD in 
various contexts and reactor configurations.

It has grown from a collection of papers that we have published over
the past 10 years.

Basic structure
---------------

aldsim provides two separate functionalities: its `core` module contains
various models implemented using nondimensional variables. These models
highlight some of the scaling laws in atomic layer deposition. 

Through its chem module, it applies the models in `core` to various idealized
models of self-limited surface kinetics.


Status
------

`aldsim` is still in development. Over the next few months it will
be expanded to incorporate a variety of models.

Quick install
-------------

Through pypi::

   pip install aldsim

Alternatively, it can be directly installed from its github repository: https://github.com/anglyan/aldsim

Acknowledgements
----------------

* Argonne Laboratory Directed Research and Development program


Copyright and license
---------------------

Copyright © 2024, UChicago Argonne, LLC

`aldsim` is distributed under the terms of BSD License. 

Argonne Patent & Intellectual Property File Number: SF-24-041


Contents
========

.. toctree::
   :maxdepth: 2

   intro
   api
   core


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
