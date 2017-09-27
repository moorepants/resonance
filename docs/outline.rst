===============
Topical Outline
===============

The course is taught over 20 two hour class periods during a quarter system of
10 weeks of instructions and 1 week for examinations. One of the 20 class
periods is reserved for a midterm examination, 2 hours are reserved for exam
reviews leaving 36 hours of in class time. The following lists the topics for
each of the class periods which correspond to the detailed headers below:

+----+----------+-----------------------+
| L# | Date     | Notebook #            |
+====+==========+=======================+
| 01 | W Sep 27 | 1, 2                  |
+----+----------+-----------------------+
| 02 | M Oct 02 | 3                     |
+----+----------+-----------------------+
| 03 | W Oct 04 | 4                     |
+----+----------+-----------------------+
| 04 | M Oct 09 | 5                     |
+----+----------+-----------------------+
| 05 | W Oct 11 | 6                     |
+----+----------+-----------------------+
| 06 | M Oct 16 | 7                     |
+----+----------+-----------------------+
| 07 | W Oct 18 | 8                     |
+----+----------+-----------------------+
| 08 | M Oct 23 | 9                     |
+----+----------+-----------------------+
| NA | T Oct 24 | Drop Date             |
+----+----------+-----------------------+
| 09 | W Oct 25 | 10                    |
+----+----------+-----------------------+
| 10 | M Oct 30 | 11                    |
+----+----------+-----------------------+
| 11 | W Nov 01 | 12                    |
+----+----------+-----------------------+
| 12 | M Nov 06 | Exam                  |
+----+----------+-----------------------+
| 13 | W Nov 08 | 13                    |
+----+----------+-----------------------+
| NA | F Nov 10 | Veterans Day Holiday  |
+----+----------+-----------------------+
| 14 | M Nov 13 | 14                    |
+----+----------+-----------------------+
| 15 | W Nov 15 | 15                    |
+----+----------+-----------------------+
| 16 | M Nov 20 | 16                    |
+----+----------+-----------------------+
| 17 | W Nov 22 | 17                    |
+----+----------+-----------------------+
| NA | R Nov 23 | Thanksgiving Holiday  |
+----+----------+-----------------------+
| NA | F Nov 24 | Thanksgiving Holiday  |
+----+----------+-----------------------+
| 18 | M Nov 27 | 18                    |
+----+----------+-----------------------+
| 19 | W Nov 29 | 19                    |
+----+----------+-----------------------+
| 20 | M Dec 04 | 20                    |
+----+----------+-----------------------+
| 21 | W Dec 06 | 21                    |
+----+----------+-----------------------+
| NA | T Dec 12 | Final Exam @ 6:00 PM  |
+----+----------+-----------------------+

Analyzing Vibrating Systems
===========================

1. Introduction to Jupyter
--------------------------

This notebook introduces students to the Jupyter notebook environment and
establishes good practices for creating computational notebooks and scientific
python programming.

After the completion of this assignment students will be able to:

- open Jupyter notebooks and operate basic functionality
- fetch assignments, complete exercises, submit work and view the graded work
- solve basic scientific python problems
- create a well formatted and fully executing notebook

2. Introduction to vibrations: Book Balancing on a Cup
------------------------------------------------------

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

3. Measuring a Bicycle Wheel's Inertia
--------------------------------------

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

4. Clock Pendulum with Air Drag Damping
---------------------------------------

This notebook introduces the third fundamental characteristic of vibration:
energy dissipation through damping. A simple pendulum model is implemented that
allows students to vary the damping parameters and visualize the three regimes
of linear damping.

After the completion of this assignment students will be able to:

- understand the concept of damped natural frequency and its relationship to
  mass/inertia, stiffness, and damping
- state the three fundamental characteristics that make a system vibrate
- compute the free response of a linear system with viscous-damping in all
  three damping regimes
- identify critically damped, underdamped, and overdamped behavior
- determine whether a linear system is over/under/critically damped given its
  dynamic properties
- understand the difference between underdamping, overdamping, and crticial
  damping

5. Clock Pendulum with Air Drag and Joint Friction
--------------------------------------------------

This notebook builds on the previous one by introducing nonlinear damping
through Coulomb friction. Students will be able to work with both a linear and
nonlinear version of the same system (pendulum) in order to compare the free
response in both cases.

After the completion of this assignment students will be able to:

- identify the function that governs the decay envelope
- compare this non-linear behavior to the linear behavior
- estimate the period of oscillation
- compute the free response of a non-linear system with viscous and coulomb
  damping

6. Vertical Vibration of a Bus Driver's Seat
--------------------------------------------

This notebook introduces external forcing of a vibratory system, where the
external force is modeled as a sinusoidal input to the bottom of a bus driver's
seat.

After the completion of this assignment students will be able to:

- excite a system with a sinusoidal input
- understand the difference in transient and steady state solutions
- use autocorrelation to determine period
- relate the frequency response to the time series
- create a frequency response plot
- define resonance and determine the parameters that cause resonance

7. Vertical vibration of a Bus Driver's Seat with a Leaf Spring
---------------------------------------------------------------

This notebook builds on the previous one by replacing the linear spring with
a realistic leaf spring.

After the completion of this assignment students will be able to:

- create a force versus displacement curve for a leaf spring
- describe the time response and frequency response of a non-linear system
- show that sinusoidal fitting does not necessarily describe non-linear
  vibration

8. Bicycle Lateral Vibration
----------------------------

This notebook introduces a simple lean and steer bicycle model as an example of
a system with multiple degrees of freedom. Coupling and modes are discussed
from a data analysis perspective.

After the completion of this assignment students will be able to:

- get a sense of the coupling of input to output through frequency response
  plots
- simulate a 2 DoF vibratory model
- identify a MDoF system and see effects of coupling through time and frequency
  domain
- determine if a general 2 DoF is stable
- sweep through input frequencies to discover modal frequencies

9. Simulating a building during an earthquake
---------------------------------------------

This notebook uses a lumped parameter multi-story building model as a
many-degree-of-freedom system with all oscillatory modes.

After the completion of this assignment students will be able to:

- examine time domain and frequency coupling with MDoF
- sweeping through frequencies to discover modal frequencies
- visualize the system's response at modal frequencies to see mode shapes

Modeling Vibrating Systems
==========================

10. Modeling the Bicycle Wheel Inertia Measurement System
---------------------------------------------------------

This notebook walks through modeling two different test rigs for determining
the vibrational characteristics of a bicycle wheel. After coming up with a
simple model the students will use the canonical linear form of the equations
of motion to derive various vibrational parameters.

After the completion of this assignment students will be able to:

- derive the equations of motion of a compound pendulum with Lagrange's method
- derive the equations of motion of a torsional pendulum with Lagrange's method
- linearize the compound pendulum equation
- put equations in canonical form
- review solutions to ODEs

11. Modeling a non-linear spring
--------------------------------

TODO : Think this out more.

After the completion of this assignment students will be able to:

- will be able to derive the nonlinear euqations of motion of a system with
  simple kinmeatics with lagrange's method

12. Modeling the car on the bumpy road
--------------------------------------

Here will will present the base excitation single degree of freedom system and
the students will derive the equations of motion. They will then explore the
displacement and force transmisiblity frequency response functions.

After the completion of this assignment students will be able to:

- derive the linear equations of motion ofa system with simple kinematics using
  lagrange's method
- create system object with custom euqations of motion an simulate the system

13. Modeling the book on a cup
------------------------------

The book balancing on the cup will be revisited. The students will derive the
equations of motion which require more complex kinematic analysis and explore
the analytical equations of motion. The stability thresholds will be determined
as well as the period from the linear model.

After the completion of this assignment students will be able to:

- derive the euqations of motion of a system with non-trivial kinematics with
  lagrange's method
- apply a linearization procedure to non-linear equations of motion
- determine the stability of a linear system analytically and verify through
  simulation

14. Balancing your car tire at the autoshop
-------------------------------------------

The mass imbalance problem will be presented through the analytical model of an
unbalance car tire. The frequency response will be derived and examined.

After the completion of this assignment students will be able to:

- derive the equations of motion fo a mass imbalance system

15. Engine cam non-sinusoidal periodic forcing
----------------------------------------------

Using an engine cam piecewise periodic function the students will learn how a
Fourier series can be used to find the solution to the differential equations
symbolicaly.

After the completion of this assignment students will be able to:

- generate a Fourier series of a periodic function
- find the analytic solution of the the mass-spring-damper system

16. Modeling a bulding during an earthquake
-------------------------------------------

We will revisit the multi-story building model and derive the equations of
motion for the system. The students will use eigenanalysis of the simple system
to discover the modes of motion and simulate the behavior.

After the completion of this assignment students will be able to:

- perform modal analysis of the system to determine its modal frequencies and
  mode shapes
- represent model using a matric equation of motion (canoncial form)
- formulate the equations of motion for a MDoF system
- use eignvalue analyssis to determine the modeshapes of a mDoF system
- plot the motion of a MDoF system (with no damping) using the analytical
  solution
- form a MDoF model corresponding to a chain of floors in a buliding

17. Bicycle Model
-----------------

The students will be given the analytical canocial form of the bicycle
equations that do not have simple damping. They will have to convert to state
space form and do a full eigenanalysis of the general form. The modes will be
examined and the nature of the bicycle motion discovered.

After the completion of this assignment students will be able to:

- convert the canonical linear form into state space form
- interpret eigenvalues and eienvectors of a general 2 DoF linear system

Designing Vibrating Systems
===========================

18. Design a Clock that Keeps Time
----------------------------------

The students will be presented with a compound pendulum model of a clock's bob
that does not keep time well due to friction and air drag. They will be tasked
with designing a system that adds in the right amount of additional energy so
that the pendulum has the desired constant period.

After the completion of this assignment students will be able to:

- develop an analytic model of a energy injection system
- simulate the motion of clock and determine its time varying period
- choose the energy injection system parameters that will cause the clock to
  work as intended

19. Isolator Selection
----------------------

The students will be presented with a model of X and asked to select and/or
design a commercially available vibration isolator that ensures the system
meets specific vibrational design criteria.

After the completion of this assignment students will be able to:

- discuss and justify trade-offs and design decisions
- model the system with additional damping provided by isolation
- select/design a vibration isolator to meet given vibration specifications
- analyze a system's motion to determine its vibrational characteristics

20. Designing a Tuned Mass Damper to Earthquake Proof a Building
----------------------------------------------------------------

Students will be presented with a single (or multi?) floor building model. They
will need to modify the model to includes a laterally actuated mass on the
roof. They will be asked to design an actuation scheme that prevents the
building from having too large of displacements or resonance while excited by a
earthquake-like vibration at its base.

After the completion of this assignment students will be able to:

- add a generic vibration absorber to a building model
- use a building model to simulate the motion of a building without damping
- choose design criteria for the building and justify decisions (with ISO
  standards)
- design an absorber that meets their design criteria
- use the frequency response function to demonstrate the effect of the
  vibration absorber

21. Designing a stable bicycle
------------------------------

The students will be presented with a 2 DoF linear model of a bicycle in
canonical form with analytical expressions for the M, C, and K  matrix entries
that are functions of the 25 bicycle parameters. The students will be asked to
discover bicycle designs that meet certain criteria through eigenanalysis and
simulation.

After the completion of this assignment students will be able to:

- determine parameters which cause the 2 DoF system to be stable/unstable
- simulate and visualize the motion of a bicycle with difference parameters
- determine and describe the influence of the physical parameters, initial
  conditions, and steering input on the dynamics of the vehicle
- design a bicycle that meets specific design criteria

22. Designing Shock Absorbtion for a Car
----------------------------------------

The students will be presented with 2D planar data generated from a "ground
truth" 3 DoF half car model. Their job will be to design a quarter car model
that behaves similarly to the ground truth model. Once they have a working
simple model, then they will design an improved shock absorber for the quarter
car model using analytic and computational methods. The instructors will then
provide the students with the ground truth model, i.e. the "real" car, and the
students will need to show that the ride quality is improved and that design
criteria is met.

After the completion of this assignment students will be able to:

- develop a simple analytic model that predicts motion provided from
  planar 2D "experimental" data
- select springs and dampers to meet given design criteria by demonstrating
  performance with the simple analytic model
- demonstrate that the designed shock absorber works well for the "real" car
- discuss why the design does or does not meet the design criteria
- reflect on their modeling and design decisions after having tested it against
  the ground truth model
