===============
Topical Outline
===============

The course is taught over 20 two hour class periods during a quarter system of
10 weeks of instructions and 1 week for examinations. One of the 20 class
periods is reserved for a midterm examination, 2 hours are reserved for exam
reviews leaving 36 hours of in class time. The following lists the topics for
each of the class periods:

01 W Sep 27
02 M Oct 02
03 W Oct 04
04 M Oct 09
05 W Oct 11
06 M Oct 16
07 W Oct 18
08 M Oct 23
NA T Oct 24 Drop Date
09 W Oct 25
10 M Oct 30
11 W Nov 01
12 M Nov 06
13 W Nov 08
NA F Nov 10 Veterans Day Holiday
14 M Nov 13
15 W Nov 15
16 M Nov 20
17 W Nov 22
NA R Nov 23 Thanksgiving Holiday
NA F Nov 24 Thanksgiving Holiday
18 M Nov 27
19 W Nov 29
20 M Dec 04
21 W Dec 06
NA T Dec 12 Final Exam @ 6:00 PM

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

19. Isolator Selection
----------------------

- discuss and justify tradeoffs and design decisions
- model the system with additional damping provided by isolation
- design a vibration isolator to meet given vibration specifications
- analyze a system's motion ot determine its vibrational characteristics

20. Designing a tuned mass damper to earthquake proof a building
----------------------------------------------------------------

- design an absorber that meets their design criteria
- choose design criteria for the buliding and justify decisions (with ISO
  standards)
- add a mass damper to the model and use the frequency repsonse function to
  demonstrate its effect
- use a buling model to simulate the motion of a building without damping

21. Designing a stable bicycle
------------------------------

- determine parameters which cause the system to be stable
- determine and descirbe the influence of the dynamic an dgoemetrical
  parameterss on stability
- explore the motion of a bicycl ewith difference dynamic parameters through
  simulation

22. Designing a shock absorbtion for a car
------------------------------------------

- use experimental data from a car to form a model
- select a damper to meet given design criteria and demonstrate this with their
  model
- validate their designed damper in the underlying (complex) car model
- discuss why their design doe or does not meet the criteria with the more
  complex model
- reflect on their modeling and design decisions after having tested it against
  the groun truth model

Implementation: give them just experimental data with a general description of
the system and inputs, have the form the model and design an absorber, only
then give them the underlying (complex) model to test their design with
