========================================================================
Resonance: Learning Mechanical Vibration Engineering Through Computation
========================================================================

.. image:: https://img.shields.io/pypi/v/resonance.svg
   :target: http://pypi.org/project/resonance

.. image:: https://readthedocs.org/projects/resonance/badge/?version=latest
   :target: http://resonance.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://travis-ci.org/moorepants/resonance.svg?branch=master
   :target: https://travis-ci.org/moorepants/resonance

Introduction
============

This repository contains the interactive learning materials designed for the
upper-level UC Davis engineering course on Mechanical Vibrations (ENG 122). The
materials are designed with these ideas in mind:

- That students can learn about mechanical vibrations engineering through
  "computational thinking" and "computational experimentation", i.e. actively
  interacting with a computer by writing code to simulate and analyze
  computational models and experimental data.
- That the computer allows students to solve vibration engineering problems
  without knowing all of the mathematical theory a priori. This means that we
  can motivate students to dig deeper into the theory and by presenting it
  posteriori when the motivation is high. The students will be introduced to
  data analysis techniques to study vibrations before analytical techniques.
- Students learn best by doing. The content is meant to used in class while the
  instructors act as a coach through the learning.
- That each lesson should have a motivated real life example that drives the
  investigation.
- Open access materials promote easy reuse, remixing, and dissemination.

The current course website can be found at:

https://moorepants.github.io/eng122/

All of the Jupyter notebooks are rendered at:

http://moorepants.github.io/resonance

Learning Objectives
===================

There are three broad learning objectives that we focus on in the course:

1. Students will be able to analyze vibrational measurement data to draw
   conclusions about the measured system's vibrational nature and describe how
   the systems behaves vibrational.
2. Students will be able to create simple mathematical and computational models
   of real vibrating systems that can be used to answer specific questions
   about the system by concisely demonstrating the vibrational phenomena.
3. Students will be able to design a mechanical structure that has desirable
   vibrational behavior.

Students that master these three core learning objectives will be well prepared
to use mechanical vibration concepts, theories, and tools to solve engineering
problems.

For a more detailed topical outline with specific per-activity learning
objectives see the `outline <outline.rst>`_.

Assessment
==========

The students will be assessed through a series of in- and out-of- class
exercises that focus on individual lesson topics, two examinations, and on an
individual open-ended vibration design project.

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

Installation
============

For users, you can create a conda environment called ``resonance`` by
downloading the ``user-environment.yml`` file and typing the following at the
command line::

   $ conda env create -f user-environment.yml

This environment can be activated with::

   $ conda activate resonance

To properly view the exercises you will need to enable the exercise2 notebook
extension::

   (resonance)$ jupyter nbextension enable exercise2/main

If you want to develop resonance, use the ``dev-environment.yml`` file::

   $ conda env create -f dev-environment.yml
   $ conda activate resonance-dev

If you don't want to use our environments, you can use pip to install
resonance::

   $ pip install resonance
