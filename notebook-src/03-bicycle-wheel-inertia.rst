======================================
3. Measuring a Bicycle Wheel's Inertia
======================================

This notebook introduces the concept of using vibratory characteristics to
estimate parameters of an existing system. It discusses how vibrations can be
measured and how these measurements might relate to parameters of interest,
such as the inertia of a bicycle wheel.

After the completion of this assignment students will be able to:

- describe different methods of measuring vibrations
- choose appropriate sensors and sensor placement
- visualize the vibrational measurements
- use curve fitting to estimate the period of oscillation
- understand the concept of natural frequency and its relationship to
  mass/inertia and stiffness
- state two of the three fundamental characteristics that govern vibration
  (mass/inertia and stiffness)
- use frequency domain techniques to characterize a system's behavior

Inertia and Vibration
=====================

One of the fundamental properties that affects how systems vibrate is the
inertia of a system. Inertia is colloquially defined as the *tendency to do
nothing or to remain unchanged*. A vibrating system is changing with respect to
time, thus we may infer that inertia will try to prevent the system from
changing. There are two specific types of inertia that we will address: mass
which is a resistance to linear motion and moments of inertia which is a
resistenace to angular motion. The moment of inertia can also be thought a
descriptor for the distribution of mass in a system.

Most people are familiar with a bicycle wheel, so we will look into the
inertial characteristics of bicycle wheels and how this relates to vibration. A
bicycle wheel is generally symmeetric about the plane defined by the circular
rim and tire in addition to being symmetric about any radial plane normal to
the wheel plane. Recalling from dynamics, this means that the momement of
inertia will have three priniciapl moments of inertia, two of which are equal.
The inertia about the axis of the wheel is the rotational inertia (resists the
wheel rolling) and the inertia about any radial axis resists motions like
turning the handlebars, etc.

The moment of inertia about the axis of a bicycle wheel can reasonably be
estimated by:

.. math:: I = \frac{m r^2}{2}

which is the theorecital value for a an infinitely thin hoop or ring. See
Wikipedia's `List of moments of inertia`_ for more.

.. _List of moments of inertia: https://en.wikipedia.org/wiki/List_of_moments_of_inertia


We demonstrated in the previous lesson that the inertia of the book affected
the frequency of oscillation.

A common way to measure inertia is to use a spring with a known linear
force/displacement ratio, the spring constant or stiffness, to resist motion.
If we attached a torsional spring along a radial line of the bicycle wheel to
the wheel and a fixed ceiling. If you twist the wheel about the axis of the
torsional spring it will start vibrating.

With a spring and a mass/inertia you can cause vibration. The spring constant
is the second fundamental property of vibration.

The video belows shows a bicycle wheel that is attached to the ceiling by a
thin steel rod. The rod's axis passing through the center of mass of the wheel.
When the bicycle wheel is given an initial angular displacement from the
equilibrium, it vibrates.

.. code:: pycon

   >>> from IPython.display import YouTubeVideo
   >>> YouTubeVideo('uuc1_RTYZLw')

A free body diagram can be sketched of this system, see below:

.. image:: fig/03/bicycle-wheel-radial-inertia-measurement-fbd.png
   :width: 400px
   :align: center

.. topic:: Exercise

   This system is slightly different than the book and cup system in the
   previous lesson. There are two fundamental properties of the system that
   make this it vibrate. What are these two things?

Data
====

During the experiment shown in the above video I recorded some important
constants about the system and measured the angular rate of the system about
the rod's axis during rotation. The data from the measurement can be loaded into a Panda's ``Series`` with Panda's ``read_csv()`` function

.. code:: pycon

   >>> import pandas as pd
   >>> gyro_reading = pd.read_csv('data/03/bicycle-wheel-radial-inertia-rate-gyro-measurement.csv', index_col='time')
   >>> gyro_reading

This can then be plotted to see what the signal looks like:

.. code:: pycon

   >>> gyro_reading.plot()

.. topic:: Exercise

   Use your period estimation function from the previous lesson to try to get an estimate of the period of this oscillation.

There is a system in the ``resonance`` package that represents the bicycle wheel and torsion rod system:

.. code:: pycon

   >>> from resonance.linear_systems import BicycleWheelRadialInertiaSystem
   >>> sys = BicycleWheelRadialInertiaSystem()

Note that the constants are not set:

.. code:: pycon

   >>> sys.constants

Here are some values for the system's geometry and mass:

- Outer radius of the bicycle wheel, :math:`r`: 0.336 m
- Mass of the bicycle wheel, :math:`m`: 1.55 kg
- Length of the torsion rod, :math:`l`: 1.05 m
- Diameter of the torsion rod, :math:`d`: 0.00635 m
- Modulus of Rigidity of steel, :math:`G`: 77 GPa

.. topic:: Exercise

   Use the above constants and simulate the system with an initial angle of 2 degrees. Be sure to use consistent units. Calculate the period of oscillation of this system and see if it is similar to the data.

Interactive Plots
=================

.. code:: pycon

   >>> import matplotlib.pyplot as plt
   >>> fig, ax = plt.subplots(1, 1)
   >>> ax.plot(rate_gyro.index, rate_gyro, label='Measured Signal')
   >>> traj = sys.free_response(10)
   >>> line = ax.plot(traj.index, traj.angle)
   >>> def plot_trajectory_comparison(radius, mass):
   ...    sys.constants['radius'] = radius
   ...    sys.constants['mass'] = mass
   ...    traj = sys.free_response(10)
   ...    line.set_data(traj.index, traj.angle)
   ...
   >>> from ipywidgets import interact
   >>> interact(plot_trajectory_comparison)

We learned in the last lesson that there is likely a relationship between the inertia of the system and the frequency of oscillation. It would be nice to plot the frequency versus the change in inertia to try and determine what the relationship is. Say we want to check a range of inertia values from X to Y, we can create those values with:

.. code::

   >>> import numpy as np
   >>> inertias = np.linspace(0.01, 0.2, num=100)

Instead of typing the simulation code out for each of the 100 inertia values, we can use a loop to iterate through each value and save some typing, e.g.:

.. code:: pycon

   >>> radii = []
   >>> for inertia in inertias:
   ...     radius = np.sqrt(2 * inertia / sys.constants['mass'])
   ...     print(radius)
   ...     radii.append(radius)
   ...

This printed the new radius at each loop iteration and also appended it to a list of radii, which is now available to use:

... code:: pycon

   >>> radii

You may also want to convert this to a NumPy array:

.. code:: pycon

   >>> radii = np.array(radii)

Note that if all you were doing was the above computation you can simply do:

.. code:: pycon

   >>> radii = np.sqrt(2 * inertias / sys.constants['mass'])

.. topic:: Exercise

   Use a loop to construct a list of frequencies for different inertia values. After you have both arrays, plot the inertias on the X axis and the frequencies on the Y axis. Is there any functional relationship that describes the relationship between the variables?


Curve Fitting

mP = 4.65+/-0.01
dP = 0.03009+/-0.00001, 0.03010+/-0.00001, 0.03012+/-0.00001
lP = 0.8355+/-0.001
period of calib rod oscillation
Tp = 0.9561002971327414

Ip = mP / 12 * (3*rp^2 + lp^2)
Ip = 0.2707614811040625
stiffness of rod
k = 4 * Ip * np.pi**2 / Tp**2
k = 11.693370530226998


k = G Jp / l
Jp = np.pi * d**4 / 32  # polar second moment of are of torsion rod
l = G * Jp / k = G * np.pi * d**4 / 32 / k
if rod is 1/4" then length of torsion rod is
1.0511042914686415 meters


IFxx = 0.0524475128396 kg m**2
IFyy = 0.0983720589324

k * T**2 / 4 / np.pi**2 = Ip
Front wheel compount pendulum length
lF = 0.2957195+/-0.00008