Learning Objectives
===================

There are three broad learning objectives that we focus on in the course:

1. Students will be able to analyze vibrational measurement data to draw
   conclusions about the measured system's vibrationnal nature and describe how
   the systems behaves vibrationally.
2. Students will be able to create simple mathematical and computational models
   of real vibrating systems that can be used to answer specific questions
   about the system by consisely demonstrating the vibrational phenomena.
3. Students will be able to design a mechanical structure that has desireable
   vibrational behavior.

Students that master these three core learning objectives will be well prepared
to use mechanical vibration concepts, theories, and tools to solve engineering
problems.

Assessment
==========

The students will be assessed through a series of in- and out-of- class
exercises that focus on individual lesson topics, two examinations, and on an
individual open-ended vibration design project.

Methodology
===========

There are a few principles to teaching:

1. We want to maximize the student "doing" both in-class and out.
2. Computational objects will be used to represent real vibrational systems and
   the students will interact with these objects to discover vibrational
   concepts as a first step.

Topical Outline
===============

The course is designed around 18 two hour in-class learning sessions, 2 two
hour exams, spread over an 11 week quarter. The class is divided into three
broad topics based on the three primary learning objectives.

Analyzing Vibrating Systems
---------------------------

After finishing this section students will be able to:

- state the three fundamental characteristics that make a system vibrate
- describe different methods of measuring vibrations
- choose appropriate sensors and a placement
- visualize the vibrational measurements
- visualize a system’s free response
- identify critically damped, underdamped, and overdamped behavior
- excite a system with different input signals
- use time domain techniques to characterize a system's behavior
- use frequency domain techniques to characterize a system's behavior
- identify a MDoF system and see effects of coupling through time and frequency domain
- interactively adjust core system parameters to affect response

The story arc:


- Start with a conversation about real things that vibrate. Show vibrating
  systems and ask for examples. Have class discussion on commonalities with the
  instructor goal of identifying mass/inertia, flexibility/stiffness, damping
  as core common attributes of vibratory systems.
- How might we measure the motion of these things? What are important
  measurements? Talk about position, velocity, and acceleration measurement
  techniques. Mass, stiffness, damping measurement?
- Introduce a 1 DoF system (block on an ice rink) (mass and damping?). Show how
  to simulate an initial value problem and have the students inspect the
  measurements (position, velocity, acceleration, .. jerk). What do you
  observe with initial velocity? (velocity decreases over time) Have them look
  at short sim, so velocity doesn't change, and then a long sim where velicity
  does change (goes to zero)
- Introduce a rotational 1 DOF system (mass and damping) and it's measurements
  (angle, angular velocity, angular acceleration)
- Now introduce a 1 DoF ?bungee jumper? system (or something without gravity).
  Explain that this sysem has an elastic element that resists the motion of
  mass. Let the students simulate this and talk about the characteristics of
  the measurements? What does position functino look like? What is the
  relationship among the three measurements? Can they determine the dreivative
  relationship?
- Show how to adjust the mass, and stiffness values of the system. Let
  students explore to discover effects of the parameters on the free response.
  Can they find out that m and k primarily affect the frequency of oscillation?
  Introduce the concept of natural frequency.
- Introduce autocorrleation to see period of signal?
- Call the signal frequency spectrum method (FFT of signal) to look at the
  frequency vs amplitude of the signal. Introduce the frequency domain view of
  the signal.
- Introduce the effects of damping on the free response: over, under, critical
- Introduce multi dof system (two bungee jumpers?), talk about dof, coordinates, etc
- Simulate initial value problem and look at signals in time and frequency
  domain.
- Now let's poke at the systems with an input. Starting with a sinusoidal
  input. What does the output look like? Should identify that it is also
  sinusoidal but diff amplitude. Try different frequencies.
- Introduce the frequency response. Use a sweep frequency input, call a
  underlying sys id method that produces the frequency response plot? Talk
  about amplitude and phase relationships.
- Explore how the mass, stifneess, and damping affect the response plot.
- Introduce resonance.
- Look at frequency response of mdof systems.

Modeling Vibrating Systems
--------------------------

After finishing this section students will be able to:

- identify the essential aspects to model (# degrees of freedom, coordinates,
  modeling assumptions)
- draw a vibration free body diagram of a real system
- write the non-linear equations of motion of a vibrating system
- write a linear equation of motion for a vibrating system
- transform a mathematical model into a computational model
- validate that the model behaves like the real system in a way that is optimal
  for answering your questions about the system

Designing Vibrating Systems
---------------------------

After finishing this section students will be able to:

- realize that mass, stiffness, damping, and geometry affect behavior in a
  complex coupled way.  Understand the tradeoffs among adjusting model
  parameters
- use common methods like isolators, mass balance, mass abosrber to control
  behavior


Notebooks
=========

Analyzing Vibrating Systems
---------------------------

1. Introduction to Jupyter
~~~~~~~~~~~~~~~~~~~~~~~~~~

- open Jupyter notebooks and operate basic functionality
- fetch assignments, complete exercises, submit work and view the graded work
- solve basic scientific python problems
- create a well formatted and fully executing notebook

2. Introduction to vibrations: Book Balancing on a Cup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- visualize a system's free response
- estimate the period of a sinusoidal vibration from a time series
- compare a computer simulation result to experimental result
- interactively adjust the book inertia to see the affect on system response
- understand the concept of natural frequency nd its relationship to
  mass/inertia

3. Measuring a Bicycle Wheel's Inertia
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- identify the fucntion that govenrs the decay envelope
- compare this non-linear behavior to the linear behavior
- estimate the period of oscillation
- compute the free response of a non-linear system with viscous and coulomb
  damping

6. Vertical Vibration of a Bus Driver's Seat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- excite a system with a sinusoidal input
- understand the difference in transient and steady state solutions
- use autocorrelation to determine period
- relate the frequence response to the time series
- create a frequency response plot
- define resonance and determine the parameters that cause resonance

7. Vertical vibration of a Bus Driver's Seat with a Leaf Spring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- create a force versus displacement curve for a leaf spring
- describe the time repsonse and frequency response of a non-linear system
- show that sinusoidal fitting does not necessarily describe non-linear
  vibration

8. Bicycle Lateral Vibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- get a sense of the couppling of input to output through frequency response
  plots
- simulate a 2 DoF vibratory model
- identify a MDoF system and see effects of couplin through time and frequency
  domain
- determine if a general 2 DoF is stable
- sweeping through frequencies to discover modal frequencies

9. Simulating a building during an earthquake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- exmaine time domain and frequency coupling with MDoF
- sweeping through frequencies to discover modal frequencies

Modeling Vibrating Systems
--------------------------

10. Modeling the Bicycle Wheel Inertia Measurement System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- derive the equations of motion of a compound pendulum with Lagrange's method
- derive the equations of motion of a torsional pendulum with Lagrange's method
- linearize the compound pendulum equation
- put equations in canoncial form
- review solutions to ODEs

11. Modeling a non-linear spring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- will be able to derive the nonlinear euqations of motion of a system with
  simple kinmeatics with lagrange's method

12. Modeling the car on the bumpy road
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- derive the linear equations of motion ofa system with simple kinematics using
  lagrange's method
- create system object with custom euqations of motion an simulate the system

13. Modeling the book on a cup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- derive the euqations of motion of a system with non-trivial kinematics with
  lagrange's method
- apply a linearization procedure to non-linear equations of motion
- determine the stability of a linear system analytically and verify through
  simulation

14. Balancing your car tire at the autoshop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- derive the equations of motion fo a mass imbalance system

15. Engine cam non-sinusoidal periodic forcing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

16. Modeling a bulding during an earthquake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- perform modal analysis of the system to determine its modal frequencies and
  mode shapes
- represent model using a matric equation of motion (canoncial form)
- formulate the equations of motion for a MDoF system
- use eignvalue analyssis to determine the modeshapes of a mDoF system
- plot the motion of a MDoF system (with no damping) using the analytical
  solution
- form a MDoF model corresponding to a chain of floors in a buliding

17. Bicycle Model
~~~~~~~~~~~~~~~~~

- convert the canonical linear form into state space form
- interpret eigenvalues and eienvectors of a general 2 DoF linear system

Designing Vibrating Systems
---------------------------

18. Design a Clock that Keeps Time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

19. Isolator Selection
~~~~~~~~~~~~~~~~~~~~~~

- discuss and justify tradeoffs and design decisions
- model the system with additional damping provided by isolation
- design a vibration isolator to meet given vibration specifications
- analyze a system's motion ot determine its vibrational characteristics

20. Designing a tuned mass damper to earthquake proof a building
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- design an absorber that meets their design criteria
- choose design criteria for the buliding and justify decisions (with ISO
  standards)
- add a mass damper to the model and use the frequency repsonse function to
  demonstrate its effect
- use a buling model to simulate the motion of a building without damping

21. Designing a stable bicycle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- determine parameters which cause the system to be stable
- determine and descirbe the influence of the dynamic an dgoemetrical
  parameterss on stability
- explore the motion of a bicycl ewith difference dynamic parameters through
  simulation

22. Designing a shock absorbtion for a car
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

18. Bumpy Road
~~~~~~~~~~~~~~

- design a car shock absorber to minimize force and displacement transmisiblity

