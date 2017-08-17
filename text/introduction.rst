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
