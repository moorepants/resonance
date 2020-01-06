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
