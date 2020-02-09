================================================================
Creating and Exercising a Custom Single Degree of Freedom System
================================================================

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`custom-sdof-system` or Jupyter notebook:
   :jupyter-download:notebook:`custom-sdof-system`.

Creating a new system
=====================

The first step is to import a "blank" ``SingleDoFLinearSystem`` and initialize
it.

.. jupyter-execute::

   from resonance.linear_systems import SingleDoFLinearSystem

   msd_sys = SingleDoFLinearSystem()

Now define the constant variables for the system. In this case, the single
degree of freedom system will be described by its mass, natural frequency, and
damping ratio.

.. jupyter-execute::

   msd_sys.constants['m'] = 1.0  # kg
   msd_sys.constants['fn'] = 1.0  # Hz
   msd_sys.constants['zeta'] = 0.1  # unitless

   msd_sys.constants

Define the coordinate and speed. The software assumes that the speed is defined
as the time derivative of the coordinate, i.e. :math:`v = \dot{x}`.

.. jupyter-execute::

   msd_sys.coordinates['x'] = 1.0  # m
   msd_sys.speeds['v'] = 0.0  # m/s

.. jupyter-execute::

   msd_sys.coordinates

.. jupyter-execute::

   msd_sys.speeds

.. jupyter-execute::

   msd_sys.states

Now that the coordinate, speed, and constants are defined the equations of
motion can be defined. For a single degree of freedom system a Python function
must be defined that uses the system's constants to compute the coefficients to
the canonical second order equation, :math:`m \dot{v} + c v + k x = 0`.

The inputs to the function are the constants. The variable names should match
those defined above on the system.

.. jupyter-execute::

   import numpy as np

   def calculate_canonical_coefficients(m, fn, zeta):
       """Returns the system's mass, damping, and stiffness coefficients given
       the system's constants."""
       wn = 2*np.pi*fn
       k = m*wn**2
       c = zeta*2*wn*m
       return m, c, k

   msd_sys.canonical_coeffs_func = calculate_canonical_coefficients

Once this function is defined and added to the system :math:`m,c,k` can be
computed using:

.. jupyter-execute::

   msd_sys.canonical_coefficients()

The period of the natural frequency can be computed with:

.. jupyter-execute::

   msd_sys.period()

All information about the system can be displayed:

.. jupyter-execute::

   msd_sys

Simulating the free response
============================

The ``free_response()`` function simulates the now fully defined system given
as an initial value problem. One or both of the coordinates and speeds must be
set to provide a free response. The following shows the response to both
:math:`x` and :math:`v` being set to some initial values.

.. jupyter-execute::

   msd_sys.coordinates['x'] = -5.0
   msd_sys.speeds['v'] = 8.0

``free_response()`` returns a Pandas ``DataFrame`` with the time values as the
index and columns for the coordinate, speed, and additionally the time
derivative of the speed (acceleration in this case). See
https://pandas.pydata.org/pandas-docs/stable/getting_started/dsintro.html for
an introduction to ``DataFrame``.

.. jupyter-execute::

   trajectories = msd_sys.free_response(5.0)
   trajectories

There are a variety of plotting methods associated with the ``DataFrame`` that
can be used to quickly plot the trajectories of the coordinate, speed, and
acceleration. See more about plotting ``DataFrames`` at
https://pandas.pydata.org/pandas-docs/stable/user_guide/visualization.html.

.. jupyter-execute::

   axes = trajectories.plot(subplots=True)

Response to change in constants
-------------------------------

This system is *parameterized* by its mass, natural frequency, and damping
ratio. It can be useful to plot the trajectories of position for different
values of :math:`\zeta` for example.

Set the initial conditions back to simply stretching the spring 1 meter:

.. jupyter-execute::

   msd_sys.coordinates['x'] = 1.0
   msd_sys.speeds['v'] = 0.0

Now change :math:`\zeta` to different values and simulate the free response to
see the different damping regimes:

Un-damped, :math:`\zeta=0`

.. jupyter-execute::

   msd_sys.constants['zeta'] = 0.0  # Unitless
   trajectories = msd_sys.free_response(5.0)
   axes = trajectories['x'].plot()

Under-damped, :math:`0<\zeta<1`

.. jupyter-execute::

   msd_sys.constants['zeta'] = 0.5  # Unitless
   trajectories = msd_sys.free_response(5.0)
   axes = trajectories['x'].plot()

Critically damped, :math:`\zeta=1`

.. jupyter-execute::

   msd_sys.constants['zeta'] = 1.0  # Unitless
   trajectories = msd_sys.free_response(5.0)
   axes = trajectories['x'].plot()

Over-damped, :math:`\zeta>1`

.. jupyter-execute::

   msd_sys.constants['zeta'] = 2.0  # Unitless
   trajectories = msd_sys.free_response(5.0)
   axes = trajectories['x'].plot()

Adding measurements
===================

It is often useful to calculate the trajectories of other quantities. Systems
in resonance allow "measurements" to be defined. These measurements are
functions of the constants, coordinates, speeds, and/or time. To create a new
measurement, create a function that returns the quantity of interest. Here a
measurement function is defined that calculates the kinetic energy
(:math:`\frac{1}{2}mv^2`) of the system and then added to the system with
variable name ``KE``.

.. jupyter-execute::

   def calculate_kinetic_energy(m, v):
      return m*v**2/2

   msd_sys.add_measurement('KE', calculate_kinetic_energy)

Once added, the measurement will be computed and added to the ``DataFrame``
containing the trajectories:

.. jupyter-execute::

   msd_sys.constants['zeta'] = 0.5  # Unitless
   trajectories = msd_sys.free_response(5.0)
   trajectories

and can be plotted like any other column:

.. jupyter-execute::

   axes = trajectories['KE'].plot()

Plotting the configuration
==========================

``resonance`` systems can plot and animate at the system's configuration. To do
so, a custom function that generates a configuration plot using matplotlib must
be defined and associated with the system. Below a plot is created to show an
orange block representing the mass and a spring attached to the block. The
``spring()`` function conveniently provides the x and y data needed to plot the
spring.

.. jupyter-execute::

   import matplotlib.pyplot as plt
   from resonance.functions import spring

   # create a new constant to describe the block's dimension, l
   msd_sys.constants['l'] = 0.2  # m

   def create_configuration_figure(x, l):

       # create a figure with one or more axes
       fig, ax = plt.subplots()

       # the `spring()` function creates the x and y data for plotting a simple
       # spring
       spring_x_data, spring_y_data = spring(0.0, x, l/2, l/2, l/8, n=3)
       lines = ax.plot(spring_x_data, spring_y_data, color='purple')
       spring_line = lines[0]

       # add a square that represents the mass
       square = plt.Rectangle((x, 0.0), width=l, height=l, color='orange')
       ax.add_patch(square)

       # add a vertical line representing the spring's attachment point
       ax.axvline(0.0, linewidth=4.0, color='black')

       # set axis limits and aspect ratio such that the entire motion will appear
       ax.set_ylim((-l/2, 3*l/2))
       ax.set_xlim((-np.abs(x) - l, np.abs(x) + l))
       ax.set_aspect('equal')

       ax.set_xlabel('$x$ [m]')
       ax.set_ylabel('$y$ [m]')

       # this function must return the figure as the first item
       # but you also may return any number of objects that you'd like to have
       # access to modify, e.g. for an animation update

       return fig, ax, spring_line, square

   # associate the function with the system
   msd_sys.config_plot_func = create_configuration_figure

Now the configuration plot can be generated with ``plot_configuration()``. This
returns the same results as the function defined above.

.. jupyter-execute::

   fig, ax, spring_line, square = msd_sys.plot_configuration()

Animating the configuration
===========================

Reset to un-damped motion and simulate again

.. jupyter-execute::

   msd_sys.constants['zeta'] = 0.1
   trajectories = msd_sys.free_response(5.0)

To animate the configuration, create a function that updates the various
matplotlib objects using any constants, coordinates, speeds, and/or the special
variable ``time``. The last input arguments to this function must be all of the
extra outputs of ``plot_configuration()`` (excluding the figure which is the
first output). The order of these must match the order of the
``plot_configuration()`` outputs.

.. jupyter-execute::

   def update_configuration(x, l, time,  # any variables you need for updating
                            ax, spring_line, square):  # returned items from plot_configuration() in same order

       ax.set_title('{:1.2f} [s]'.format(time))

       xs, ys = spring(0.0, x, l/2, l/2, l/8, n=3)
       spring_line.set_data(xs, ys)

       square.set_xy((x, 0.0))

   msd_sys.config_plot_update_func = update_configuration

Now that the update function is associated, ``animate_configuration()`` will
create the animation. Here the frames-per-second are set to an explicit value.

.. jupyter-execute::

   animation = msd_sys.animate_configuration(fps=30)

If using the notebook interactively with ``%matplotlib widget`` set, the
animation above will play. But ``animate_configuration()`` returns a matplotlib
``FuncAnimation`` object which has other options that allow the generation of
different formats, see
https://matplotlib.org/api/_as_gen/matplotlib.animation.FuncAnimation.html for
options. One option is to create a Javascript/HTML versions that displays
nicely in the notebook with different play options:

.. jupyter-execute::

   from IPython.display import HTML

   HTML(animation.to_jshtml(fps=30))

Response to sinusoidal forcing
==============================

The response to a sinusoidal forcing input, i.e.:

.. math::

   m\dot{v} + cv + kx = F_o \sin(\omega t)

can be simulated with ``sinusoidal_forcing_response()``. This works the same as
``free_response`` except it requires a forcing amplitude and frequency.

.. jupyter-execute::

   msd_sys.coordinates['x'] = 0.0  # m
   msd_sys.speeds['v'] = 0.0  # m/s

   Fo = 10.0
   omega = 2*np.pi*3.0  # rad/s

   forced_trajectory = msd_sys.sinusoidal_forcing_response(Fo, omega, 5.0)

Note that there is now a ``forcing_function`` column. This is the applied
forcing function.

.. jupyter-execute::

   forced_trajectory

The trajectories can be plotted and animated as above:

.. jupyter-execute::

   axes = forced_trajectory.plot(subplots=True)

.. jupyter-execute::

   fps = 30
   animation = msd_sys.animate_configuration(fps=fps)

.. jupyter-execute::

   HTML(animation.to_jshtml(fps=fps))

Frequency response
==================

The frequency response to sinusoidal forcing at different frequencies can be
plotted with ``frequency_response_plot()`` for a specific forcing amplitude.

.. jupyter-execute::

   axes = msd_sys.frequency_response_plot(Fo)

Response to periodic forcing
============================

Any periodic forcing function can be applied given the Fourier series
coefficients of the approximating function. The following function calculates
the Fourier series coefficients for a "sawtooth" shaped periodic input.

.. jupyter-execute::

   def sawtooth_fourier_coeffs(A, N):
       """
       A : sawtooth amplitude, Newtons
       T : sawtooth period, seconds
       N : number of Fourier series terms
       """
       n = np.arange(1, N+1)
       an = A*(8*(-1)**n - 8) / 2 / np.pi**2 / n**2
       return 0, an, np.zeros_like(an)

   a0, an, bn = sawtooth_fourier_coeffs(Fo, 20)

These coefficients can be provided to ``periodic_forcing_response()`` to
simulate the response:

.. jupyter-execute::

   wb = 2*np.pi*3.0  # rad/s

   trajectory = msd_sys.periodic_forcing_response(a0, an, bn, wb, 5.0)
   trajectory

.. jupyter-execute::

   axes = trajectory.plot(subplots=True)

.. jupyter-execute::

   fps = 30
   animation = msd_sys.animate_configuration(fps=fps)

.. jupyter-execute::

   HTML(animation.to_jshtml(fps=fps))
