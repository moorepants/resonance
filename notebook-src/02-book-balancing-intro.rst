======================================================
2. Introduction to vibrations: Book Balancing on a Cup
======================================================

This notebook introduces a single degree of freedom vibratory system in which a
textbook balances on a cylindrical cup. The system is implemented as a model
that students can interact with in order to visualize its free response and
compare to the demonstration in the classroom.

After the completion of this assignment students will be able to:

- visualize a system's free response
- estimate the period of a sinusoidal vibration from a time series
- compare a computer simulation result to experimental result
- interactively adjust the book inertia to see the affect on system response
- understand the concept of natural frequency nd its relationship to
  mass/inertia

What are vibrations?

Here we will study a simple vibratory system. If you set a cylindrical cup on
its side and place a book on top of it the book will oscillate if initialized
at a small angle non-horizontal angle or it is perturbed slightly. Note that it
oscillates about an equilibrium point, the rest angle defines this equilibrium
point. Vibration is defined as an oscialation about an equilibrium.

During this class we will examine and explore many different vibrating systems,
such as this simple book and cup system. We will have some live demos as we are
showing now but in general we will work with computational systems to
experiment and learn about a system. For example, here is a system that
represents the book on a cup::

   >>> from resonance.linearsystems import BookOnCupSystem

The line above loads the command that can create systems that behave like the
live demo at teh fron of the room. To create a system::

   >>> sys = BookOnCupSystem()

Systems have different attributes, for example this system has geometry, such
as the books height, width, and length and the cup's radius. The book also has
a mass and in this case assumed to be uniformly dense. You can set the values
of these attributes, as such::

   >>> from resonance.units import centimeter, kilogram
   >>> sys.parameters.height.value = 1 * centimeter
   >>> sys.parameters.width.value = 6 * centimeter
   >>> sys.parameters.length.value = 10 * centimeter
   >>> sys.parameters.radius.value = 3 * centimeter
   >>> sys.parameters.mass.value = 1 * kilogram

TODO : Should a parameter be a subclass of SymPy symbol with values, units, and
a long form name?

::

   >>> sys.parameters.radius
   Parameter(r, 'radius', value=0.03 * meter)

Above we import some units from resonance and specify the value of the system's
parameters. All systems will have different sets of parameters. This system may
have more parameters too, for example what if you were on the moon? Maybe the
acceelration due to gravity could be changed or what if the book to cup
connection was very slipperly, maybe the coefficient of friction would be a
specific parameter. Parameters, as we'e defined here are constant with respect
to time.

There are other types of parameters too. A very important type of parameter are
those that varying with time. What are the time varying parameters of this
system?

Two types states and measurements.

TODO : Access states via dictionary key or attribute?
TODO : Should a time varying parameter be able to store an array of values associated with an array of times?

.. codeblock:: pycon

   >>> sys.state['book angle']
   TimeVaryingParameter(theta, 'book angle', value=0 * radian, t=0 * second)
   >>> sys.state.book_angular_rate
   TimeVaryingParameter(omega, 'book_angular_rate', value=0 * radian / second, t=0 * second)
   # TODO : Should units be separate than the value? Units don't work well with numpy/scipy.
   >>> sys.state.book_angular_rate.value
   0 * radian / second
   >>> sys.state.book_angular_rate.units
   radian / second

There are also measurements. These measurements are, in general a function of
the states and the parameters. For example, maybe you really want to know the
height of the top edge of the book as as function of time. A measurement can be
defined as such::

   >>> from sympy import cos
   >>> theta = sys.state['book angle']
   >>> h = sys.par['height']
   >>> w = sys.par['width]
   >>> end_height = h * cos(theta) + sin(theta) * w / 2
   >>> sys.measurements['end_height'] = end_height
   >>> sys.measurements['end_height'].value()  # computes the value based on current value of parameters and states
   >>> 1 * centimeter

If you change the book angle you'll get a different measurement::

   >>> from resonance.units import degree
   >>> sys.state['book angle'].value = 1 * degree
   >>> sys.measurements['end_height'].value()

TODO : Explain the difference in value as a attribute and value as a method.

Initial Value Problem
=====================

Now that we have a system with defined constant parameters we can make it move,
in our case vibrate. There are two ways to create motion: apply perturbing
forces to the system or set the state values to an initial value other than the
equilibrium state. We will do the later here. We can set the initial angle to 1
degree and then simulate the system::

   >>> sys.state['book angle'] = 1 * degree
   >>> trajectory = sys.simulate(t0=0 * second, tf=5 * second)
   >>> trajectory
             book angle  book angular rate  end height
   time
   0.000000    0.000000           1.000000    0.000000
   0.555556    0.527415           0.849608    0.527415
   1.111111    0.896192           0.443666    0.896192
   1.666667    0.995408          -0.095724    0.995408
   2.222222    0.795220          -0.606321    0.795220
   2.777778    0.355842          -0.934546    0.355842
   3.333333   -0.190568          -0.981674   -0.190568
   3.888889   -0.679658          -0.733529   -0.679658
   4.444444   -0.964317          -0.264750   -0.964317
   5.000000   -0.958924           0.283662   -0.958924

TODO : Is a pandas data frame a good thing to return here? Not sure if we
should stick to common python objects or have custom ones. Allen may have a
custom time series type object for his class. One issue with having time as the
index is that we'd like to represent as a float but indices shouldn't be
represented as a float.

TODO : Should the trajectory be stored on the system object so that other
methods can access and use the most recent simulation? For example, what if I
want to call the FFT method on one column? Or do we let them use the FFT numpy
functions manually?

We can now plot these variables, one at a time::

   >>> trajectory['book angle'].plot()

altogether::

   >>> trajectory.plot()

or in subplots::

   >>> trajectory.plot(subplots=True)

Exercise
--------

Have them simulate and plot different initial conditions. They should be able
to identify what functions govern the motion.


Time Series Analaysis
=====================

Find the period of oscillation
Autocorrelation?
