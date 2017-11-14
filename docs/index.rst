.. resonance documentation master file, created by
   sphinx-quickstart on Tue Sep 26 12:16:37 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to resonance's documentation!
=====================================

Resonance is a companion software library to an interactive textbook written for
an upper level undergraduate introduction to mechanical vibrations course.

The library is specifically designed for learning engineering principles
through computational thinking and computational experimentation and thus we
have guidelines on its design to facilitate this. Those guidelines are roughly:

- Don’t teach programming for the sake of teaching programming. Show them how
  to solve problems and introduce programming along the way to solve those
  problems.
- Hide the fine details of programming and only use simple constructs, so that
  learning vibrations is highlighted instead of programming.
- Hide the simulation details (linear/nonlinear ODE solutions).
- Centered around the “System” object. Systems represent real things: a car, a
  bridge, a bicycle, an airplane wing.
- Students can use and construct systems.
- Students only create functions, no need to understand classes and their
  construction.
- Easy visualizations (time history plots and animations of systems)
- Extra informative and lots of error messages (try to predict student
  mistakes)

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   outline.rst
   api.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
