=========================================================================================
2. What Are Mechanical Vibrations?: Investigating a Book Oscillating on a Cylindrical Cup
=========================================================================================

Introduction
============

This notebook introduces a single degree of freedom vibratory system in which a
textbook balances on a cylindrical cup. The system is implemented as a model
that students can interact with in order to visualize its free response and
compare the computer simulation to the demonstration in the classroom.

After the completion of this assignment students will be able to:

- view and set parameters on a system
- visualize a system's free response
- estimate the period of a sinusoidal vibration from a time series
- compare a computer simulation result to experimental result
- interactively adjust the book inertia to see the affect on system response
- understand the concept of natural frequency nd its relationship to
  mass/inertia

What are vibrations?
====================

Here we will study a simple vibratory system. A vibrating mechanical system is
typically defined as a collection of rigid and flexible objects that interact
in a closed envelope. If you set a cylindrical cup on its side and place a book
on top of it, the book will oscillate if initially displaced at a small
non-horizontal angle. Note that it oscillates about a horizontal position, i.e.
the rest angular position is called an equilibrium point, equilibrium state, or
equilibrium configuration. Vibration is formally defined as an oscillation
about an equilibrium.

During this class we will examine and explore many different vibratory systems,
such as this simple book and cup system. We will have some live demos, as we
are showing now, but in general we will work with computational representations
of systems to experiment and learn about the nature of vibration. For example,
here is a system that represents the book on a cup::

   >>> from resonance.linear_systems import BookOnCupSystem

The line above loads the command that can create systems that behave like the
live demo at the front of the room. To create a system, execute the following
cell by pressing the shift and enter key simultaneously::

   >>> sys = BookOnCupSystem()

Systems have different parameters, for example this system has geometry, such
as the book's height, width, and length and the cup's radius. The book also has
a mass and in this case assumed to be uniformly dense. You can view all of the
parameters, which are stored in a Python dictionary by typing::

   >>> sys.parameters
   {'height': 0.029, 'length': 0.238, 'radius': 0.042, 'mass': 1.058}

A Python dictionary maps keys to values. For example the key ``'height'`` is
associated with a value ``0.029``. An individual parameter value can be
accessed by using square brackets and the key like so::

   >>> sys.parameters['radius']
   0.042
   >>> type(sys.parameters['radius'])
   float

These should be Python floating point numbers. You can set the values of these
attributes, as such::

   >>> sys.parameters['height'] = 1.0  # cm
   >>> sys.parameters['width'] = 6.0  # cm
   >>> sys.parameters['length'] = 10.0  # cm
   >>> sys.parameters['radius'] = 3.0  # cm
   >>> sys.parameters['mass'] = 1.0  # kg

Note that you will be responsible for ensuring that the units are consistent
and that all angles should be in radians. Load the system again to get back the
default parameters.

.. code-block:: pycon

   >>> sys = BookOnCupSystem()

All systems will have different sets of parameters. This system may have more
parameters too, for example what if you were on the moon? Maybe the
acceleration due to gravity could be changed or what if the book to cup
connection was very slippery, maybe the coefficient of friction would be a
specific parameter. Parameters, as we've defined here, are constant with
respect to time.

Time Varying Parameters
=======================

There are other types of parameters too. A very important type of parameter are
those that vary with time.

Exercise
--------

What are the time varying parameters of this system?


There are are an infinite number of time varying parameters. Most systems are
described by a uniquely simple set of time varying parameters, often called
generalized coordinates. These coordinates define the configuration of the
system. In our case the vertical and horizontal location of the book's mass
center could uniquely describe the configuration of the system (if the book
can't slip on the cup). But a better choice would be to use the single time
varying angle of the books surface relative to horizontal to define the
configuration. The angle of the book is a generalized coordinate because no
fewer number of parameters can possible be used to describe the configuration.
This number of generalized coordinates corresponds to the number of degrees of
freedom of a system. The degrees of freedom are XXX. The non-slipping book on a
cup has 1 degree of freedom which is described by the single generalized
coordinate, the book's angle. The system's generalized coordinates can be
accessed as such::

   >>> sys.coordinates
   {'book_angle': 0.0}
   >>> sys.coordinates['book_angle']
   0.0

Another type of time varying parameter that can be extracted from systems are
*measurement parameters*. For example, maybe you are interested in the vertical
and horizontal location of the book's center of mass instead of the generalized
coordinate. These two measurement parameters are a function of the book angle
and the system's geometry. We will use Python functions to define this
relationship. Given the value of a generalized coordinate and the values of the
system's parameters, we can define a function that computes the measurement
parameter. For example:

.. code-block:: pycon

   >>> import numpy as np
   >>> def compute_y_mass_location():
   ...     # in the sys class this function will be wrapped and all of the
   ...     # parameters, coordinates, and measurements will be injected into the
   ...     # namespace just above the function so the students can just write
   ...     # these functions with the correct variables and not worry about
   ...     # unpacking arrays
   ...     return ((radius + height / 2) * np.cos(book_angle) + radius *
   ...              book_angle * np.sin(book_angle) + radius)
   ...
   >>> sys.add_measurement('mass_center_height', compute_vertical_mass_location)

.. code-block:: pycon

   >>> def bottom_left_y(radius, height, length, book_angle):
   ...     r = radius
   ...     h = height
   ...     l = length
   ...     theta = book_angle
   ...     return r + r * np.cos(theta) + (r * theta + l / 2) * np.sin(theta)
   ...
   >>> sys.add_meas('bottom_left_y', bottom_left_y)

.. code-block:: pycon

   >>> def bottom_left_x(radius, height, length, book_angle):
   ...     r = radius
   ...     h = height
   ...     l = length
   ...     theta = book_angle
   ...     return r * np.sin(theta) - (r * theta + l / 2) * np.cos(theta)
   ...
   >>> sys.add_meas('bottom_left_x', bottom_left_x) 

TODO : Explain a Python function.

If you change the book angle you'll get a different measurement:

.. code-block:: pycon

   >>> sys.coordinates['book_angle'] = np.deg2rad(1)
   >>> # calls __getitem__ of a Measurements class and compute the right value using the supplied function
   >>> sys.measurements['mass_center_height']
   5.8

Exercise
--------

Create a measurement for the horizontal position of the center of mass of the
book and call it ``mass_center_y``.

Initial Value Problem
=====================

Now that we have a system with defined constant parameters we can make it move,
in our case vibrate. There are two ways to create motion: apply perturbing
forces to the system or set the state values to an initial value other than the
equilibrium state. We will do the later here. We can set the initial angle to 1
degree and then simulate the system::

   >>> sys.coordinates['book_angle'] = np.deg2rad(1)
   >>> trajectories = sys.simulate(t0=0, tf=5)

This creates what is called a DataFrame. DataFrames are defined in the Pandas
Python package and are essentially 2D tables with labels for each column and an
index for each row. In our case the index is the time value and the columns are
the coordinates and the measurements::

   >>> type(trajectories)
   DataFrame
   >>> trajectories
             book_angle  mass_center_height  bottom_left_x  bottom_left_y
   Time [s]
   0.000000    0.017453            0.098504      -0.118982       0.086083
   0.016722    0.017322            0.098504      -0.118982       0.086067
   0.033445    0.016929            0.098504      -0.118983       0.086021
   0.050167    0.016282            0.098504      -0.118984       0.085943
   0.066890    0.015389            0.098503      -0.118986       0.085836
   0.083612    0.014264            0.098503      -0.118988       0.085702
   0.100334    0.012925            0.098502      -0.118990       0.085541
   0.117057    0.011390            0.098502      -0.118992       0.085358
   0.133779    0.009684            0.098501      -0.118994       0.085154
   0.150502    0.007832            0.098501      -0.118996       0.084933
   0.167224    0.005862            0.098500      -0.118998       0.084698
   0.183946    0.003804            0.098500      -0.118999       0.084453
   0.200669    0.001689            0.098500      -0.119000       0.084201
   0.217391   -0.000452            0.098500      -0.119000       0.083946
   0.234114   -0.002587            0.098500      -0.119000       0.083692
   0.250836   -0.004682            0.098500      -0.118999       0.083443
   0.267559   -0.006706            0.098501      -0.118997       0.083203
   0.284281   -0.008630            0.098501      -0.118996       0.082975
   0.301003   -0.010424            0.098501      -0.118994       0.082762
   0.317726   -0.012060            0.098502      -0.118991       0.082568
   0.334448   -0.013515            0.098503      -0.118989       0.082396
   0.351171   -0.014766            0.098503      -0.118987       0.082247
   0.367893   -0.015795            0.098503      -0.118985       0.082126
   0.384615   -0.016586            0.098504      -0.118984       0.082032
   0.401338   -0.017127            0.098504      -0.118983       0.081968
   0.418060   -0.017409            0.098504      -0.118982       0.081935
   0.434783   -0.017430            0.098504      -0.118982       0.081932
   0.451505   -0.017188            0.098504      -0.118982       0.081961
   0.468227   -0.016687            0.098504      -0.118983       0.082020
   0.484950   -0.015934            0.098503      -0.118985       0.082109
   ...              ...                 ...            ...            ...
   4.515050   -0.003055            0.098500      -0.118999       0.083637
   4.531773   -0.005137            0.098500      -0.118998       0.083389
   4.548495   -0.007142            0.098501      -0.118997       0.083151
   4.565217   -0.009039            0.098501      -0.118995       0.082926
   4.581940   -0.010801            0.098502      -0.118993       0.082717
   4.598662   -0.012399            0.098502      -0.118991       0.082528
   4.615385   -0.013810            0.098503      -0.118989       0.082361
   4.632107   -0.015014            0.098503      -0.118987       0.082218
   4.648829   -0.015991            0.098504      -0.118985       0.082103
   4.665552   -0.016727            0.098504      -0.118983       0.082015
   4.682274   -0.017212            0.098504      -0.118982       0.081958
   4.698997   -0.017437            0.098504      -0.118982       0.081932
   4.715719   -0.017399            0.098504      -0.118982       0.081936
   4.732441   -0.017099            0.098504      -0.118983       0.081971
   4.749164   -0.016541            0.098504      -0.118984       0.082037
   4.765886   -0.015735            0.098503      -0.118985       0.082133
   4.782609   -0.014691            0.098503      -0.118987       0.082256
   4.799331   -0.013425            0.098502      -0.118989       0.082406
   4.816054   -0.011958            0.098502      -0.118992       0.082580
   4.832776   -0.010310            0.098501      -0.118994       0.082775
   4.849498   -0.008507            0.098501      -0.118996       0.082989
   4.866221   -0.006576            0.098501      -0.118997       0.083218
   4.882943   -0.004546            0.098500      -0.118999       0.083459
   4.899666   -0.002447            0.098500      -0.119000       0.083709
   4.916388   -0.000312            0.098500      -0.119000       0.083963
   4.933110    0.001829            0.098500      -0.119000       0.084218
   4.949833    0.003941            0.098500      -0.118999       0.084469
   4.966555    0.005995            0.098500      -0.118998       0.084714
   4.983278    0.007958            0.098501      -0.118996       0.084948
   5.000000    0.009801            0.098501      -0.118994       0.085168

   [300 rows x 4 columns]

The result of the last simulation is always stored on the system for later use.
Data frames have a ``head()`` function that shows just the first lines::

   >>> sys.results.head()
            book_angle  mass_center_height  bottom_left_x  bottom_left_y
   Time [s]
   0.000000    0.017453            0.098504      -0.118982       0.086083
   0.016722    0.017322            0.098504      -0.118982       0.086067
   0.033445    0.016929            0.098504      -0.118983       0.086021
   0.050167    0.016282            0.098504      -0.118984       0.085943
   0.066890    0.015389            0.098503      -0.118986       0.085836

::

   >>> %matplotlib inline

We can now plot these variables, one at a time::

::

   >>> trajectories['book_angle'].plot()

altogether::

   >>> trajectories.plot()

or in subplots::

   >>> trajectories.plot(subplots=True)

Maybe you want to use degrees instead, just make a new column::

   >>> trajectories['book_angle_deg'] = np.rad2deg(trajectories['book_angle'])
   >>> trajectories['book_angle_deg'].plot()

Exercise
--------

Simulate the system with different initial conditions and parameter values.

- Does the simulation always work, if not what doesn't work? *Hint: try a tall
  stack of books, can you find a stack height that is significant?*
- Are there any mathematical functions that could be used describe the change
  in the book angle?

Animate The Motion
==================

Plotting the coordinates and measurements as a function with respect to time is
a very useful way to visualize a system's motion, but it is often quite helpful
to animate a pictorial diagram of the system for easier visualization of the
motion. matplotlib has

::

   >>> import matplotlib.pyplot as plt
   >>> from matplotlib.patches import Circle, Rectangle
   >>> def fig_setup(time, radius, length, height, bottom_left_x, bottom_left_y)):
   ...     fig, ax = plt.subplots(1, 1)
   ...     ax.set_xlim((-0.15, 0.15))
   ...     ax.set_ylim((0.0, 0.2))
   ...     ax.set_xlabel('x [m]')
   ...     ax.set_ylabel('y [m]')
   ...     ax.set_aspect('equal')
   ...
   ...     circ = Circle((0.0, radius), radius=radius)
   ...
   ...     rect = Rectangle((bottom_left_x, bottom_left_y),
   ...                      length, height,
   ...                      angle=-np.rad2deg(book_angle),
   ...                      color='black')
   ...
   ...     ax.add_patch(circ)
   ...     ax.add_patch(rect)
   ...
   ...     text = ax.text(-0.125, 0.025, 'Time = {:0.3f} s'.format(time))
   ...
   ...     # return the figure first followed by any objects that change during the animation
   ...     return fig, circ, rect, text
   >>>
   >>> def animate(time, book_angle, bottom_left_x, bottom_left_y):
   ...
   ...     text.set_text('Time = {:0.3f} s'.format(time))
   ...
   ...     rect.set_xy((bottom_left_x, bottom_left_y))
   ...
   ...     # TODO : This should be a public set_angle method.
   ...     rect._angle = -np.rad2deg(book_angle)
   >>> sys.configuration_plot_function = figure_setup
   >>> sys.configuration_plot_update_function = animate
   >>> sys.plot_configuration()
   >>> sys.animate_configuration()

Exercise
--------

Using different initial conditions and parameters, compare the animation with
the time series plots.

Exercies
--------

Using the ``ax.set_title()`` function, make the title display the time value of
time so that it updates with the correct time during each animation frame.

Time Series Analysis
====================

From the above plots you can see that the oscillation is periodic and for most
cases sinusoidal. Using your program, create a function that calculates the
period of the non-linear model to three significant figures of the 11
oscillations when the initial book angle is X degrees. Compare the period
predicted by the system to the period measured in class.

Hint: Look for sign changes with np.sign(), use boolean indexing to extract
important times, and finally np.diff() and np.mean() can be useful for finding
the delta times and averaging. Note that np.diff() returns one fewer item in
the array it operates on.

::

   def find_period(t, theta):
       """Computes the period of oscillation based on the trajectory of theta.

       Parameters
       ==========
       t : array_like, shape(n,)
           An array of monotonically increasing time values.
       theta : array_like, shape(n,)
           An array of values for theta at each time in ``t``.

       Returns
       =======
       T : float
           An estimate of the period of oscillation.

       """

       peak_idxs = np.diff(np.sign(theta)) < 0
       peak_idxs = np.hstack((peak_idxs, False))
       T = np.diff(t[peak_idxs]).mean()

       return T
