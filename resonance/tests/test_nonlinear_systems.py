import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from ..nonlinear_systems import SimplePendulum

sys = SimplePendulum()

sys.constants['mass'] = 1.0  # kg
sys.constants['length'] = 1.0  # m
sys.constants['acc_due_to_grav'] = 9.81  # m/s**2

sys.coordinates['angle'] = np.deg2rad(15.0)  # rad
sys.speeds['angle_vel'] = 0.0  # rad/s


def equations_of_motion(state, time, length, acc_due_to_grav):

    l = length
    g = acc_due_to_grav

    theta = state[0]
    omega = state[1]

    thetad = omega
    omegad = -g / l * np.sin(theta)

    return [thetad, omegad]

# function that computes the right hand side of the equations of motion in
# first order form
sys.equations_of_motion = equations_of_motion


def height(angle, length):
    """The Y coordinate points in the opposite of gravity, i.e. up. The X
    coordinate points to the right."""
    return -length * np.cos(angle)

sys.add_measurement('height', height)


def sway(angle, length):
    """The X coordinate points to the right."""
    return length * np.sin(angle)

sys.add_measurement('sway', sway)


def potential_energy(length, height, mass, acc_due_to_grav):
    return mass * acc_due_to_grav * (length + height)

sys.add_measurement('potential_energy', potential_energy)


def kinetic_energy(mass, length, angle_vel):
    return mass * (length * angle_vel)**2 / 2.0

sys.add_measurement('kinetic_energy', kinetic_energy)


def total_energy(kinetic_energy, potential_energy):
    return kinetic_energy + potential_energy

sys.add_measurement('total_energy', total_energy)

trajectories = sys.free_response(5.0)


def plot_config(time, mass, length, sway, height, potential_energy):
    fig, axes = plt.subplots(1, 2)

    circle_radius = mass / 2.0 / 10.0

    axes[0].set_xlim((-length - circle_radius, length + circle_radius))
    axes[0].set_ylim((-length - circle_radius, length + circle_radius))
    axes[0].set_aspect('equal')

    rod_lines = axes[0].plot([0, sway], [0, height])

    circle = Circle((sway, height), radius=circle_radius)

    axes[0].add_patch(circle)

    axes[1].set_ylim((0, 0.5))
    axes[1].set_xlim((0, 5))
    pe_lines = axes[1].plot([time], [potential_energy], marker='o',
                            markersize=10)

    return fig, circle, rod_lines, pe_lines

sys.config_plot_func = plot_config


def update_plot(sway, height, circle, time, potential_energy, rod_lines,
                pe_lines):
    circle.center = sway, height
    rod_lines[0].set_data([0, sway], [0, height])
    pe_lines[0].set_data([time], [potential_energy])

sys.config_plot_update_func = update_plot
