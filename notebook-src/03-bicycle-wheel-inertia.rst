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

The inertia of a bicycle wheel can be estimated by:

.. math:: I = m r^2

A common way to measure inertia is to use a spring with a known linear
force/displacement ratio, the spring constant or stiffness, to resist motion.
If we attached a torsional spring along a radial line of the bicycle wheel to
the wheel and a fixed ceiling. If you twist the wheel about the axis of the
torsional spring it will start vibrating.

With a spring and a mass/inertia you can cause vibration. The spring constant
is the second fundamental property of vibration.


