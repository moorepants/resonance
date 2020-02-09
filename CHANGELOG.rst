0.21.0
======

- Switched to MIT license.
- Corrected the equation (and output) of ``periodic_forcing_response()``:
  https://github.com/moorepants/resonance/pull/206
- Added documentation example of creating a single degree of freedom linear
  system: https://github.com/moorepants/resonance/pull/205
- Added better __repr__ formats for the custom ordered dictionaries:
  https://github.com/moorepants/resonance/pull/204
- Added a __repr__ for System: https://github.com/moorepants/resonance/pull/203
- Added several notebooks from the Winter 2020 offering of ENG 122.
- Added jupyter-sphinx as a dependency to build the documentation.

0.20.0
======

- Bumped up lowest supported dependencies (approximately those in Ubuntu 18.04
  LTS). Python 3.6-3.8 now supported and code adjusted for deprecated features.
  Note that, due to conda having trouble installing the 18.04 version set, the
  tests are run with slightly newer NumPy, Pandas, and SciPy but should work
  with the 18.04 set.
- Updated doctr details for automatic builds of the notebooks on Travis CI.
- Added centered rectangle and spring functions for assisting in plotting and
  animations.
- Added a spring to the quarter car animation.
- Added the 2D "Bicycle Model" of a car.
