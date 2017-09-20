=============================================================================================
2. Introducing Mechanical Vibrations By Investigating a Book Oscillating on a Cylindrical Cup
=============================================================================================

Introduction
============

This notebook introduces a single degree of freedom vibratory system of a
textbook balancing on a cylindrical cup. The system is implemented as a model
that students can interact with in order to visualize its free response and
compare the computer simulation to a demonstration in the classroom.

Learning Objectives
-------------------

After the completion of this assignment students will be able to:

- load a specific system
- view and set parameters on a system
- create a visualization that depicts the system's configuration
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
lying flat on top of it, the book will oscillate if initially displaced at a
small non-horizontal angle. Note that it oscillates about a horizontal
position. This position is called an *equilibrium point*, *equilibrium state*,
or *equilibrium configuration* which is a natural position the system comes to
when there is no motion. Vibration is formally defined as an *oscillation about
an equilibrium*.

During this class we will examine and explore many different vibratory systems,
such as this simple book on a cup system. We will have some live demos, as we
are showing now, but in general we will work with computational representations
of systems to experiment and learn about the nature of vibration. For example,
here is a system that represents the book on a cup.  To execute the following
Jupyter code cell press the shift and enter key simultaneously:

.. code:: pycon

   >>> from resonance.linear_systems import BookOnCupSystem

The line above loads the command, ``BookOnCupSystem``, that can create systems
that behave like the live demo at the front of the room. To create a system and
assign it to a variable called ``sys``, execute the following cell:

.. code:: pycon

   >>> sys = BookOnCupSystem()

This system has many variables and functions associated with it. You can see
what the system has and can do by using the Python dot notation. Type ``sys.``
and press the tab key to see a list of all the variables and functions that the
system applies.

Constant Parameters
===================

Systems have different parameters, for example this system has geometry, such
as the book's height and length and the cup's radius. The book also has a mass
and, in this case, an underlying assumption is that the book is uniformly
dense. Note that these parameters do not change with time, i.e. they are
constant with respect to time. You can view all of the constant parameters,
which are stored in a Python dictionary by typing:

.. code:: pycon

   >>> sys.parameters
   {'height': 0.029, 'length': 0.238, 'radius': 0.042, 'mass': 1.058}

A Python dictionary maps keys, the variable name, to values, the default
numerical value. For example the key ``'height'`` is associated with a value
``0.029``. An individual parameter value can be accessed by using square
brackets:

.. code:: pycon

   >>> sys.parameters['radius']
   0.042

You can set the values of these attributes as such:

.. code:: pycon

   >>> sys.parameters['height'] = 1.0  # cm
   >>> sys.parameters['width'] = 6.0  # cm
   >>> sys.parameters['length'] = 10.0  # cm
   >>> sys.parameters['radius'] = 3.0  # cm
   >>> sys.parameters['mass'] = 1.0  # kg

*Note that you will be responsible for ensuring that the units are consistent
and that all angles should be in radians.* Overwrite the ``sys`` variable by
loading the system again to get back the default parameters.

.. code:: pycon

   >>> sys = BookOnCupSystem()
   >>> sys.parameters
   {'height': 0.029, 'length': 0.238, 'radius': 0.042, 'mass': 1.058}

All systems will have different sets of constant parameters. This system could
also have more parameters, for example what if you were on the moon? Maybe the
acceleration due to gravity could be changed. Or what if the book to cup
connection was very slippery? Maybe the coefficient of friction would be a
specific parameter. It is important to note that the system parameters, as
we've defined here, are constant with respect to time.

Time Varying Parameters
=======================

There are other types of parameters too. Another very important type of
parameter are those that vary with time.

Exercise
--------

What are the time varying parameters of this system?


There are are an infinite number of time varying parameters. Most systems are
described by a uniquely simple set of time varying parameters, often called
*generalized coordinates*. These coordinates define the configuration of the
system. In our case the vertical and horizontal location of the book's mass
center could uniquely describe the configuration of the system (if the book
can't slip on the cup). But a better choice would be to use the single time
varying angle of the books surface relative to horizontal to define the
configuration. The angle of the book is thus a generalized coordinate because
no fewer number of time varying parameters can possibly be used to describe the
configuration.  For simple systems, the number of generalized coordinates
corresponds to the number of *degrees of freedom* of a system. The degrees of
freedom are the number of independent parameters that define the configuration.
The non-slipping book on a cup has 1 degree of freedom which is described by
the single generalized coordinate, the book's angle. The system's generalized
coordinates can be accessed as such:

.. code:: pycon

   >>> sys.coordinates
   {'book_angle': 0.0}
   >>> sys.coordinates['book_angle']
   0.0

.. _degrees of freedom: https://en.wikipedia.org/wiki/Degrees_of_freedom_(mechanics)

Another type of time varying parameter that can be extracted from systems are
non-generalized coordinates or other combinations of constant parameters. We
will call these *measurement parameters*. For example, maybe you are interested
in the vertical and horizontal location of the book's lower left corner. If you
had an appropriate dynamics distance sensor you could measure this location as
it changes with time. Or you can define how these two measurement parameters
are a function of the book angle and the system's geometry. You can add
measurement parameters to a system by defining a Python functions that compute
the geometric relationship. Given the value of a generalized coordinate and the
values of the system's constant parameters, we can define a function that
computes the measurement parameter. These functions need to be to accept
scalars and arrays, so import NumPy to make use of array aware functions like
``cos()`` and ``sin()``:

.. code:: pycon

   >>> import numpy as np

The height of the bottom left corner of the book relative to the origin is
defined by this mathematical function:

.. math::

   P_y = r + r \cos{\theta} + (r \theta + l / 2) \sin{\theta}

and the Python function that implements this would look like:

.. code:: pycon

   >>> def bottom_left_y(radius, height, length, book_angle):
   ...     r = radius
   ...     h = height
   ...     l = length
   ...     theta = book_angle
   ...     return r + r * np.cos(theta) + (r * theta + l / 2) * np.sin(theta)
   ...

Note that the variable names in the function signature must be exactly as you
have defined them in the parameters. This function can be used independently
like so:

.. code:: pycon

   >>> angles = np.deg2rad(np.arange(10))
   >>> bottom_left_y(0.042, 0.029, 0.238, angles)

And you can add this measurement function to the system, so that it
automatically calculates the measurement for you with:

.. code:: pycon

   >>> sys.add_measurement('bottom_left_y', bottom_left_y)

Similarly, you can add the horizontal position:

.. code:: pycon

   >>> def bottom_left_x(radius, height, length, book_angle):
   ...     r = radius
   ...     h = height
   ...     l = length
   ...     theta = book_angle
   ...     return r * np.sin(theta) - (r * theta + l / 2) * np.cos(theta)
   ...
   >>> sys.add_measurement('bottom_left_x', bottom_left_x)

Now, if you change the book angle you'll get a different measurement:

.. code:: pycon

   >>> sys.coordinates['book_angle'] = np.deg2rad(1)
   >>> sys.measurements['bottom_left_y']
   5.8

Notice that all angles must be in radians and that we use the NumPy function
``rad2deg()`` to do this.

Visualizing the Configuration
=============================

It is often very helpful to visualize a system's configuration. In this case we
need a two dimensional drawing similar to the diagram above. The package
``matplotlib`` provides "patch" objects that represent a circle and a rectangle
that should be able to make a reasonable diagram. First import these functions:

.. code:: pycon

   >>> import matplotlib.pyplot as plt
   >>> from matplotlib.patches import Circle, Rectangle

And for nice interactive plotting in the Jupyter notebook, use this command to
turn it on:

.. code:: pycon

   >>> %matplotlib inline

Just as we did with the measurement functions, you can create a function that
generates the figure using the system's various parameters.

.. code:: pycon

   >>> def create_plot(radius, length, height, book_angle, bottom_left_x, bottom_left_y):
   ...     # create a blank figure and set basic settings on the axis
   ...     fig, ax = plt.subplots(1, 1)
   ...     ax.set_xlim((-0.15, 0.15))
   ...     ax.set_ylim((0.0, 0.2))
   ...     ax.set_xlabel('x [m]')
   ...     ax.set_ylabel('y [m]')
   ...     ax.set_aspect('equal')
   ...
   ...     # circles are created by supplying an (x, y) pair and the radius
   ...     circ = Circle((0.0, radius), radius=radius)
   ...
   ...     # rectangles are created by supplying the (x, y) pair locating the
   ...     # bottom left corner, the width, the height, and the to rotation
   ...     # angle. notice that the rotation angle is defined in the opposite
   ...     # direction as we have and it is supposed to be in degrees not radians
   ...     rect = Rectangle((bottom_left_x, bottom_left_y),
   ...                      length, height,
   ...                      angle=-np.rad2deg(book_angle),
   ...                      color='black')
   ...
   ...     ax.add_patch(circ)
   ...     ax.add_patch(rect)
   ...
   ...     return fig
   ...

A system can have a single configuration plot function and it is set like so:

.. code:: pycon

   >>> sys.config_plot_func = create_plot

Now have a look at the configuration plot by calling ``plot_configuration()``:

.. code:: pycon

   >>> sys.plot_configuration();

Exercise
--------

Change the system's parameters (constant parameters and the book angle) and see
how the plot reflects these changes.

Free Response
=============

Now that we have a system with defined constant parameters we can make it
vibrate. There are two ways to create this motion: apply perturbing forces to
the system or set the coordinate to an initial angle other than the equilibrium
angle. We will do the later here. The resulting motion is called the *free
response* of the system, meaning that no external forces are causing the
motion. To simulate the free response of the system, some values of time are
needed. In this case an initial time value and a final time value are passed
into the ``free_response()`` function. First, set the initial angle of the book
and then call ``free_repsonse)``, storing the result in a variable named
``trajectories``:

.. code:: pycon

   >>> sys.coordinates['book_angle'] = np.deg2rad(1)
   >>> trajectories = sys.free_response(0, 5)

This creates what is called a data frame. Data frames are defined in the Pandas
Python package and are one of the most common Python data types. They are
essentially a 2D table with labels for each column and an index for each row.
In our case the index is the time value and the columns are the values of the
coordinates and the measurements at that particular time:

.. code:: pycon

   >>> type(trajectories)
   pandas.core.frame.DataFrame
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

The data frames have useful plotting functions associated with them, so it is
rather easy to plot the various coordinates and measurements versus time:

.. code:: pycon

   >>> trajectories.plot();

It is often a little clearer if each column is in a subplots, especially if
they have different units, as these do:

.. code:: pycon

   >>> trajectories.plot(subplots=True);

A single column can be accessed and plotted too:

.. code:: pycon

   >>> trajectories['book_angle'].plot();

Maybe you want to use degrees for the book angle instead, just make a new
column:

.. code:: pycon

   >>> trajectories['book_angle_deg'] = np.rad2deg(trajectories['book_angle'])
   >>> trajectories['book_angle_deg'].plot();

Exercise
--------

Create the free response of the system with different initial coordinate values
and parameter values.

- Does the simulation always work, if not what doesn't work? *Hint: try a tall
  stack of books, can you find a stack height that is significant?*
- What mathematical function can be used describe the change in the book angle?

Animate The Motion
==================

Now that we we have a time varying response, we can animate the configure
figure to visualize how the system moves. There is one minor change that needs
to be made to the configuration plot function first. We need to make sure that
it returns any of the objects that change with time.

.. code:: pycon

   >>> def create_plot(radius, length, height, book_angle, bottom_left_x, bottom_left_y):
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
   ...     # make sure to return the rectangle, which moves at each time step!
   ...     return fig, rect
   ...
   >>> sys.config_plot_func = create_plot

Now, an animation update function can be created which updates the bottom left
corner's x and y coordinate at each time step. The last argument in the
function signature must be the object that changes.

.. code:: pycon

   >>> def update_frame(book_angle, bottom_left_x, bottom_left_y, rect):
   ...     rect.set_xy((bottom_left_x, bottom_left_y))
   ...     rect._angle = -np.rad2deg(book_angle)
   ...

Lastly, add this function to the system:

.. code:: pycon

   >>> sys.config_plot_update_func = update_frame

The visualization can now be animated with:

.. code:: pycon

   >>> %matplotlib notebook
   >>> sys.animate_configuration()

Exercie
-------

There is a special variable ``time`` that can be specified in the plot setup
and update functions. Add this variable to the function signatures and create
some text on the plot that displays the current time using:

.. code:: python

   # text = ax.text(-0.125, 0.025, 'Time = {:0.3f} s'.format(time))
   # text.set_text('Time = {:0.3f} s'.format(time))

Time Series Analysis
====================

From the above plots you can see that the oscillation is periodic and
sinusoidal. Using your program, create a function that calculates the period of
the non-linear model to three significant figures of the 11 oscillations when
the initial book angle is X degrees. Compare the period predicted by the system
to the period measured in class.

*Hint: Look for sign changes with np.sign(), use boolean indexing to extract
important times, and finally np.diff() and np.mean() can be useful for finding
the delta times and averaging. Note that np.diff() returns one fewer item in
the array it operates on.*

.. code:: pycon

   >>> def find_period(t, theta):
   ...     """Computes the period of oscillation based on the trajectory of theta.
   ...
   ...     Parameters
   ...     ==========
   ...     t : array_like, shape(n,)
   ...         An array of monotonically increasing time values.
   ...     theta : array_like, shape(n,)
   ...         An array of values for theta at each time in ``t``.
   ...
   ...     Returns
   ...     =======
   ...     T : float
   ...         An estimate of the period of oscillation.
   ...
   ...     """
   ...
   ...     peak_idxs = np.diff(np.sign(theta)) < 0
   ...     peak_idxs = np.hstack((peak_idxs, False))
   ...     T = np.diff(t[peak_idxs]).mean()
   ...
   ...     return T
   ...
   >>> find_period(trajectories.index, trajectories.book_angle)

Exercise
--------

Plot the period versus change in mass, length, and radius. Is there anything
interesting about these plots?
