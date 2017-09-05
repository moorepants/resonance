===============
Topical Outline
===============

The course is taught over 20 two hour class periods during a quarter system of
10 weeks of instructions and 1 week for examinations. One of the 20 class
periods is reserved for a midterm examination, 2 hours are reserved for exam
reviews leaving 36 hours of in class time. The following lists the topics for
each of the class periods:

| 01 W Sep 27
| 02 M Oct 02
| 03 W Oct 04
| 04 M Oct 09
| 05 W Oct 11
| 06 M Oct 16
| 07 W Oct 18
| 08 M Oct 23
| NA T Oct 24 Drop Date
| 09 W Oct 25
| 10 M Oct 30
| 11 W Nov 01
| 12 M Nov 06
| 13 W Nov 08
| NA F Nov 10 Veterans Day Holiday
| 14 M Nov 13
| 15 W Nov 15
| 16 M Nov 20
| 17 W Nov 22
| NA R Nov 23 Thanksgiving Holiday
| NA F Nov 24 Thanksgiving Holiday
| 18 M Nov 27
| 19 W Nov 29
| 20 M Dec 04
| 21 W Dec 06
| NA T Dec 12 Final Exam @ 6:00 PM

Analyzing Vibrating Systems
===========================

1. Introduction to Jupyter
--------------------------

- open Jupyter notebooks and operate basic functionality
- fetch assignments, complete exercises, submit work and view the graded work
- solve basic scientific python problems
- create a well formatted and fully executing notebook

2. Introduction to vibrations: Book Balancing on a Cup
------------------------------------------------------

- visualize a system's free response
- estimate the period of a sinusoidal vibration from a time series
- compare a computer simulation result to experimental result
- interactively adjust the book inertia to see the affect on system response
- understand the concept of natural frequency nd its relationship to
  mass/inertia

3. Measuring a Bicycle Wheel's Inertia
--------------------------------------

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

- understand the concept of damped natural frequency and its relationship to
  mass/inertia, stiffness, and damping
- state the three fundamental characteristics that make a system vibrate
- compute the free response of a linear system with viscous-damping in all
  three damping regimes
- identify critically damped, underdamped, and overdamped behavior
- determine whether a system is over/under/critically damped given its dynamic
  properties
- understnad the difference between underdamping, overdamping, and crticial
  damping

5. Clock Pendulum with Air Drag and Joint Friction
--------------------------------------------------

- identify the fucntion that govenrs the decay envelope
- compare this non-linear behavior to the linear behavior
- estimate the period of oscillation
- compute the free response of a non-linear system with viscous and coulomb
  damping

6. Vertical Vibration of a Bus Driver's Seat
--------------------------------------------

- excite a system with a sinusoidal input
- understand the difference in transient and steady state solutions
- use autocorrelation to determine period
- relate the frequence response to the time series
- create a frequency response plot
- define resonance and determine the parameters that cause resonance

7. Vertical vibration of a Bus Driver's Seat with a Leaf Spring
---------------------------------------------------------------

- create a force versus displacement curve for a leaf spring
- describe the time repsonse and frequency response of a non-linear system
- show that sinusoidal fitting does not necessarily describe non-linear
  vibration

8. Bicycle Lateral Vibration
----------------------------

- get a sense of the couppling of input to output through frequency response
  plots
- simulate a 2 DoF vibratory model
- identify a MDoF system and see effects of couplin through time and frequency
  domain
- determine if a general 2 DoF is stable
- sweeping through frequencies to discover modal frequencies

9. Simulating a building during an earthquake
---------------------------------------------

- exmaine time domain and frequency coupling with MDoF
- sweeping through frequencies to discover modal frequencies

Modeling Vibrating Systems
==========================

10. Modeling the Bicycle Wheel Inertia Measurement System
---------------------------------------------------------

- derive the equations of motion of a compound pendulum with Lagrange's method
- derive the equations of motion of a torsional pendulum with Lagrange's method
- linearize the compound pendulum equation
- put equations in canoncial form
- review solutions to ODEs

11. Modeling a non-linear spring
--------------------------------

- will be able to derive the nonlinear euqations of motion of a system with
  simple kinmeatics with lagrange's method

12. Modeling the car on the bumpy road
--------------------------------------

- derive the linear equations of motion ofa system with simple kinematics using
  lagrange's method
- create system object with custom euqations of motion an simulate the system

13. Modeling the book on a cup
------------------------------

- derive the euqations of motion of a system with non-trivial kinematics with
  lagrange's method
- apply a linearization procedure to non-linear equations of motion
- determine the stability of a linear system analytically and verify through
  simulation

14. Balancing your car tire at the autoshop
-------------------------------------------

- derive the equations of motion fo a mass imbalance system

15. Engine cam non-sinusoidal periodic forcing
----------------------------------------------

16. Modeling a bulding during an earthquake
-------------------------------------------

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
