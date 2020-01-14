import numpy as np
import sympy as sm
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from resonance.nonlinear_systems import SingleDoFNonLinearSystem


def print_eq(lhs, rhs):
    """Returns a SymPy equality of the provided left and right hand sides.

    Parameters
    ==========
    lhs : string
        A valid LaTeX string.
    rhs : sympifiable
        A SymPy expression.

    """
    return sm.Eq(sm.Symbol(lhs), rhs)


def eom_in_first_order_form(T, U, q, u, t):
    """Returns the equations of motion of a system.abs

    Parameters
    ==========
    T
    U
    q
    u
    t
    """
    L = T - U
    eom = L.diff(q.diff()).diff(t) - L.diff(q)
    mass = eom.expand().coeff(q.diff(t, 2))
    kin = u
    dyn = -eom.xreplace({q.diff(t, 2): 0}).xreplace({q.diff(): u}) / mass
    dyn = sm.simplify(dyn)
    return sm.Eq(sm.Matrix([q.diff(), u.diff()]), sm.Matrix([kin, dyn]))



class BookCupSystem(SingleDoFNonLinearSystem):

    def __init__(self):

        super(BookCupSystem, self).__init__()

        self.constants['d'] = 0.029  # m
        self.constants['l'] = 0.238  # m
        self.constants['r'] = 0.042  # m
        self.constants['m'] = 1.058  # kg
        self.constants['g'] = 9.81  # m/s**2
        self.coordinates['theta'] = 0.0  # rad
        self.speeds['omega'] = 0.0  # rad/s
        
        def rhs(theta, omega, d, l, r, g):
            """Returns the derivatives of the state variables.
            
            Parameters
            ==========
            theta : float
                Angle of the book in radians.
            omega : float
                Angular rate of the book in radians per second.
            d : float
                Book thickness in meters.
            l : float
                Book width in meters.
            r : float
                Cup radius in meters.
            g : float
                Acceleration due to gravity in meters per squared seconds.
                
            Returns
            =======
            thetadot : float
                The angular rate of the book in radians per second.
            omegadot : float
                The angular acceleration of the book in radians per second.
            
            """
            thetadot = omega
            omegadot = (6 * d * g * np.sin(theta) - 12 * g * r * theta * np.cos(theta) -
                12 * r**2 * omega**2 * theta) / (4 * d**2 + l**2 + 12 * r**2 * theta**2)
            return thetadot, omegadot

        self.diff_eq_func = rhs

        def bottom_left_x(r, l, theta):
            return r * np.sin(theta) - (r * theta + l / 2) * np.cos(theta)

        self.add_measurement('bottom_left_x', bottom_left_x)

        def bottom_left_y(r, l, theta):
            return r + r * np.cos(theta) + (r * theta + l / 2) * np.sin(theta)

        self.add_measurement('bottom_left_y', bottom_left_y)

        def create_plot(time, r, l, d, theta, bottom_left_x, bottom_left_y):
            fig, ax = plt.subplots(1, 1)
            width = max(l, r * 2)
            ax.set_xlim((-width / 2 - width / 10, width / 2 + width / 10))
            ax.set_ylim((0.0, 2 * r + 2 * d))
            ax.set_xlabel('x [m]')
            ax.set_ylabel('y [m]')
            ax.set_aspect('equal')

            circ = pat.Circle((0.0, r), radius=r)

            rect = pat.Rectangle((bottom_left_x, bottom_left_y),
                                 l, d,
                                 angle=-np.rad2deg(theta),
                                 color='black')

            ax.add_patch(circ)
            ax.add_patch(rect)

            title = ax.set_title('Time: {:.2f} [s]'.format(time))

            # make sure to return the rectangle, which moves at each time step!
            return fig, rect, title

        self.config_plot_func = create_plot

        def update_frame(time, theta, bottom_left_x, bottom_left_y, rect,
                         title):
            title.set_text('Time: {:.2f} [s]'.format(time))
            rect.set_xy((bottom_left_x, bottom_left_y))
            rect.angle = -np.rad2deg(theta)
            rect._angle = -np.rad2deg(theta)

        self.config_plot_update_func = update_frame
