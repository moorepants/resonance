========================================================================
Resonance: Learning Mechanical Vibration Engineering Through Computation
========================================================================

**NOTE : This material is very much in draft form. Comments welcome in the
issues.**

Introduction
============

This repository contains the interactive textbook designed for the upper-level
UC Davis engineering course on Mechanical Vibrations (ENG 122). The textbook is
designed with these ideas in mind:

- That students can learn about mechanical vibrations engineering through
  "computational thinking", i.e. actively interacting with a computer and
  writing code to simulate and analyze computational models and experimental
  data.
- That the computer allows students to solve vibration engineering problems
  without knowing all of the mathematical theory a priori. This means that we
  can motivate students to dig deeper into the theory and by presenting it
  posteriori when the motivation is high.
- Students learn best by doing. The content is meant to used in class while the
  instructors act as a coach through the learning.
- Open access materials promote easy reuse, remixing, and dissemination.

The course website can be found at https://moorepants.github.io/eng122/

Authors
=======

- Jason K. Moore, Faculty, Mechanical and Aerospace Engineering Department,
  University of California, Davis
- Kenneth Lyons, Graduate Student, Mechanical and Aerospace Engineering
  Department, University of California, Davis

License
=======

The contents of this repository are licensed under the CC-BY 4.0 license.

Acknowledgements
================

Much of this work has been made possible through the Undergraduate
Instructional Innovation Program funds provided by the Association of American
Universities (AAU) and Google which is administered by UC Davis's Center for
Educational Effectiveness.

This work is also made possible by the broad open source software stack that
underpins the Scientific Python Ecosystem, in particular: Jupyter, NumPy,
SymPy, SciPy, and matplotlib.

Learning Objectives
===================

There are three broad learning objectives that we focus on in the course:

1. Students will be able to analyze vibrational measurement data to draw
   conclusions about the measured system's vibrationnal nature and describe how
   the systems behave vibrationally.
2. Students will be able to create simple mathematical and computational models
   of real vibrating systems that can be used to answer specific questions
   about the system by consisely demonstrate the vibrational phenomena.
3. Students will be able to design a mechanical structure that has desireable
   vibrational behavior.

These core 

After completing the course, students will be able to:

- model a physical vibrating system and understand it's essential motion characteristics
- numerically integrate non-linear systems
- understand how the eigenvalue problem can model a vibrational system
- design a mechanical system with vibration in mind
- measure vibrations and interpret the meaning of the data
- to analyze a vibrational system using computation

Outline
=======

The course is taught over 20 two hour class periods during a quarter system of
10 weeks of instructions and 1 week for examinations. One of the 20 class
periods is reserved for a midterm examination, 2 hours are reserved for exam
reviews leaving 36 hours of in class time. The following lists the topics for
each of the class periods:

1 Introduction [100 mins]
-------------------------

The goal here is to lead the students through solving a vibration engineering
problem while simultaneously introducing the foundation Python commands they
will need to build on. This shouldn't be overwhelming with too many new Python
things but should ease them in but show them that they can solve something
useful the first day.

- Students will be able to use the core Python computing commands in the class.
- Students Will be able solve an introductory vibration engineering problem.

TODO : What problem to open with?

2 Modeling Vibrating Systems [50 min]
-------------------------------------

The key thing here is that we want students to be able to look at real physical
objects and visualize what the essential motion is. They should be able to
sketch out free body diagrams that indicate:

- an appropriate number of degrees of freedom
- appropriate generalized coordinate definitions
- appropriate lumped elements that describe the rigid bodies and lumped
  elements (springs, dampers) and external loads acting on the system

I'd also like them to come away with an appreciation of why it may be important
to try to create the simplest model that is capable of explaining the phenomena
of interest.

3 Formulating Equations of Motion [150 min]
-------------------------------------------

- recognize that the relationship between mass/inertia, acceleration, and the
  loads acting on a system are second order ordinary differential equations
- be able to express the linear and angular velocity magnitude expressions of
  important points and rotating reference frames
- be able to write the system's kinetic energy in terms of generalized
  coordinates
- be able to write the system's potential energy in terms of the generalized
  coordinates
- be able to form the Lagrangian
- be able write the Lagrange's equation of the first kind and evaluate it
- add non-conservative forces with Rayleigh's principle
- recognize advantages over a Newton-Euler formulation
- be able to convert between first and second order form
- be able to convert from canonical form to state space
- single and multi dof
- linearizing eoms

https://en.wikipedia.org/wiki/Lagrangian_mechanics

4 Free harmonic motion with and without viscous damping
-------------------------------------------------------

- determine the natural frequency
- find the solution to the ODE
- compare different ways to express ODE solution
- write the SDoF system in terms of nat. freq and damping ratio
- underdamped, critically damped, overdamped

5. Estimating system parameters from vibrations (live experiment where we give them data)
6. Forced harmonic motion with and without viscous damping
7. Non-linear vibration (Coulomb) + simulation of non-linear systems
8. Impulse response (heaviside)
9. Stability: book balance
10. Base excitation: car on bumpy road
11. Mass imbalance
12. Arbitrary forcing (convolution integral)
13. Arbitrary periodic forcing (Fourier series)
14. Modal analysis of decomposable systems: building
15. Modal analysis of non-decomposable: bicycle modeshapes
16. Isolator design
17. Vibration absorbers

Other:

- Stiffness
- Equivalency in stiffness, damping, mass, etc
- Free response to two dof
- Transform methods
- Response random inputs
- Analogy to electrical circuits or other energy domains
- More in-depth non-linear vibratory systems
- Relationship to FEA of structures
- Beams and membranes: continuous systems (Euler’s beam equation)
