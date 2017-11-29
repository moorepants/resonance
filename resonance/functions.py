import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.animation as animation


def estimate_period(time, signal):
    """Computes the period of oscillation based on the given periodic signal.

    Parameters
    ==========
    time : array_like, shape(n,)
        An array of monotonically increasing time values.
    signal : array_like, shape(n,)
        An array of values for the periodic signal at each time in ``t``.

    Returns
    =======
    period : float
        An estimate of the period of oscillation.

    """
    peak_idxs = np.diff(np.sign(signal)) < 0
    peak_idxs = np.hstack((peak_idxs, False))
    period = np.diff(time[peak_idxs]).mean()

    return period


def benchmark_par_to_canonical(p):
    """
    Returns the canonical matrices of the Whipple bicycle model linearized
    about the upright constant velocity configuration. It uses the parameter
    definitions from [Meijaard2007]_.
    Parameters
    ----------
    p : dictionary
        A dictionary of the benchmark bicycle parameters. Make sure your units
        are correct, best to ue the benchmark paper's units!
    Returns
    -------
    M : ndarray, shape(2,2)
        The mass matrix.
    C1 : ndarray, shape(2,2)
        The damping like matrix that is proportional to the speed, v.
    K0 : ndarray, shape(2,2)
        The stiffness matrix proportional to gravity, g.
    K2 : ndarray, shape(2,2)
        The stiffness matrix proportional to the speed squared, v**2.
    """
    mT = p['mR'] + p['mB'] + p['mH'] + p['mF']
    xT = (p['xB'] * p['mB'] + p['xH'] * p['mH'] + p['w'] * p['mF']) / mT
    zT = (-p['rR'] * p['mR'] + p['zB'] * p['mB'] +
          p['zH'] * p['mH'] - p['rF'] * p['mF']) / mT

    ITxx = (p['IRxx'] + p['IBxx'] + p['IHxx'] + p['IFxx'] + p['mR'] *
            p['rR']**2 + p['mB'] * p['zB']**2 + p['mH'] * p['zH']**2 + p['mF'] *
            p['rF']**2)
    ITxz = (p['IBxz'] + p['IHxz'] - p['mB'] * p['xB'] * p['zB'] -
            p['mH'] * p['xH'] * p['zH'] + p['mF'] * p['w'] * p['rF'])
    p['IRzz'] = p['IRxx']
    p['IFzz'] = p['IFxx']
    ITzz = (p['IRzz'] + p['IBzz'] + p['IHzz'] + p['IFzz'] +
            p['mB'] * p['xB']**2 + p['mH'] * p['xH']**2 + p['mF'] * p['w']**2)

    mA = p['mH'] + p['mF']
    xA = (p['xH'] * p['mH'] + p['w'] * p['mF']) / mA
    zA = (p['zH'] * p['mH'] - p['rF'] * p['mF']) / mA

    IAxx = (p['IHxx'] + p['IFxx'] + p['mH'] * (p['zH'] - zA)**2 +
            p['mF'] * (p['rF'] + zA)**2)
    IAxz = (p['IHxz'] - p['mH'] * (p['xH'] - xA) * (p['zH'] - zA) + p['mF'] *
            (p['w'] - xA) * (p['rF'] + zA))
    IAzz = (p['IHzz'] + p['IFzz'] + p['mH'] * (p['xH'] - xA)**2 + p['mF'] *
            (p['w'] - xA)**2)
    uA = (xA - p['w'] - p['c']) * np.cos(p['lam']) - zA * np.sin(p['lam'])
    IAll = (mA * uA**2 + IAxx * np.sin(p['lam'])**2 +
            2 * IAxz * np.sin(p['lam']) * np.cos(p['lam']) +
            IAzz * np.cos(p['lam'])**2)
    IAlx = (-mA * uA * zA + IAxx * np.sin(p['lam']) + IAxz *
            np.cos(p['lam']))
    IAlz = (mA * uA * xA + IAxz * np.sin(p['lam']) + IAzz *
            np.cos(p['lam']))

    mu = p['c'] / p['w'] * np.cos(p['lam'])

    SR = p['IRyy'] / p['rR']
    SF = p['IFyy'] / p['rF']
    ST = SR + SF
    SA = mA * uA + mu * mT * xT

    Mpp = ITxx
    Mpd = IAlx + mu * ITxz
    Mdp = Mpd
    Mdd = IAll + 2 * mu * IAlz + mu**2 * ITzz
    M = np.array([[Mpp, Mpd], [Mdp, Mdd]])

    K0pp = mT * zT  # this value only reports to 13 digit precision it seems?
    K0pd = -SA
    K0dp = K0pd
    K0dd = -SA * np.sin(p['lam'])
    K0 = np.array([[K0pp, K0pd], [K0dp, K0dd]])

    K2pp = 0.
    K2pd = (ST - mT * zT) / p['w'] * np.cos(p['lam'])
    K2dp = 0.
    K2dd = (SA + SF * np.sin(p['lam'])) / p['w'] * np.cos(p['lam'])
    K2 = np.array([[K2pp, K2pd], [K2dp, K2dd]])

    C1pp = 0.
    C1pd = (mu * ST + SF * np.cos(p['lam']) + ITxz / p['w'] *
            np.cos(p['lam']) - mu*mT*zT)
    C1dp = -(mu * ST + SF * np.cos(p['lam']))
    C1dd = (IAlz / p['w'] * np.cos(p['lam']) + mu * (SA +
            ITzz / p['w'] * np.cos(p['lam'])))
    C1 = np.array([[C1pp, C1pd], [C1dp, C1dd]])

    return M, C1, K0, K2


class Phasor(object):
    """Phasor that can be advanced in time with rotation and growth rates.

    Parameters
    ----------
    init : complex
        Initial phasor in rectangular form (Re + jIm)
    frequency : float, optional
        Rotation rate in rad/s.
    growth_rate : float, optional
        Exponential growth rate (decay if < 0).

    Attributes
    ----------
    t : float
        Current time.
    re : float
        Current real component of the phasor.
    im : float
        Current imaginary component of the phasor.
    radius : float
        Current radius of the phasor.
    angle : float
        Current angle of the phasor.
    trace_t : list
        History of time values (since most recent `clear()`).
    trace_re : list
        History of real component values (since most recent `clear()`).
    trace_im : list
        History of imaginary component values (since most recent `clear()`).
    """

    def __init__(self, init, frequency=0, growth_rate=0):
        self.init = init
        self.frequency = frequency
        self.growth_rate = growth_rate

        self.trace_t = []
        self.trace_re = []
        self.trace_im = []

        self._init()

    @classmethod
    def from_eig(cls, eigvec_component, eigval):
        """Creates a phasor from an eigenvalue/eigenvector component pair.

        Parameters
        ----------
        eigvec_component : complex
            A single eigenvector component representing the phasor's initial
            real/imaginary parts.
        eigval : complex
            The eigenvector, which specifies the phasor's growth rate (real
            part) and rotational frequency (imaginary part).
        """
        return cls(eigvec_component, np.imag(eigval), np.real(eigval))

    def advance(self, dt):
        """Advance the phasor by a time step dt."""
        self.t += dt
        self.radius *= np.exp(self.growth_rate * dt)
        self.angle += self.frequency * dt

        self._update_rect()

        self.trace_t.append(self.t)
        self.trace_re.append(self.re)
        self.trace_im.append(self.im)

    def clear(self):
        """Clear trajectories."""
        self._init()
        del self.trace_t[:]
        del self.trace_re[:]
        del self.trace_im[:]

    def _init(self):
        """Initialize parameters."""
        self.t = 0
        self.radius = np.abs(self.init)
        self.angle = np.angle(self.init)
        self._update_rect()

    def _update_rect(self):
        """Update rectangular components."""
        vec = self.radius * np.exp(1j * self.angle)
        self.re = np.real(vec)
        self.im = np.imag(vec)


class PhasorAnimation(animation.TimedAnimation):
    """Animation for demonstrating rotating phasors.

    Two axes are set up. On top, there is an s-plane to show the real and
    imaginary components of the phasors. The current phasor "vector" is shown
    with a thick line, the current endpoint of the vector is shown with a
    circle, thin lines show the projection of the real part of the phasor down
    to the bottom of the plane, and the time history of the endpoint of the
    vectors are shown.

    On bottom, the phasors' real components are plotted in time. The plot is
    rotated so that time is positive downward, and the x axes of the s-plane
    and the time plots are lined up. The current value is shown with a circle,
    thin lines show the projection from the top of the plot to the current
    value, and the time history is plotted.

    Parameters
    ----------
    fig : Figure
        matplotlib Figure object on which to animate.
    t : array
        Array of time values at which to plot. Even time spacing is assumed.
    phasors : list
        List of Phasor objects to advance and plot.
    re_range : tuple, optional
        Limits of the real axis.
    im_range : tuple, optional
        Limits of the imaginary axis.
    repeat : bool, optional
        Specifies whether or not to repeat the animation once it finishes.
    repeat_delay : float, optional
        Amount of time to wait before repeating the animation in milliseconds.
    time_stretch : float, optional
        Multiplicative factor of the plotting interval. Increasing
        `time_stretch` effectively makes the animation slower without affecting
        the time units.
    blit : bool, optional
        Specifies whether or not to use blitting.
    """

    def __init__(self, fig, t, phasors, re_range=(-1, 1), im_range=(-1, 1),
                 repeat=True, repeat_delay=0, time_stretch=1, blit=True):
        self.t = t
        self.dt = t[1] - t[0]
        self.phasors = phasors
        self.re_range = re_range
        self.im_range = im_range

        gs = GridSpec(2, 1, height_ratios=[1, 2])

        # s-plane plot of phasors
        ax_s = fig.add_subplot(gs[0])
        ax_s.set_ylabel('Im')
        ax_s.set_xlim(*re_range)
        ax_s.set_ylim(*im_range)
        ax_s.set_xticklabels(ax_s.get_xticklabels(), visible=False)
        ax_s.grid()

        # time plot of the real part of the phasors
        ax_t = fig.add_subplot(gs[1])
        ax_t.set_xlabel('Re')
        ax_t.set_ylabel('t')
        ax_t.set_xlim(*re_range)
        ax_t.set_ylim(self.t[0], self.t[-1]+self.dt)
        ax_t.invert_yaxis()
        ax_t.grid()

        fig.subplots_adjust(hspace=0)
        fig.tight_layout()

        # vectors in the Re/Im axis from origin
        self.vec_lines = []
        # dots at the end of the vectors in the Re/Im axis
        self.vec_dots = []
        # lines streaming from the endpoints of the vectors to the axis base
        self.vec_connector_lines = []
        # trace of vectors in Re/Im axis
        self.vec_trace_lines = []
        # dot showing current x(t) value
        self.time_dots = []
        # lines streaming from the current x(t) value to the top of the axis
        self.time_connector_lines = []
        # trace of x(t) values
        self.time_trace_lines = []

        color_vals = np.linspace(0, 1, 10)
        for i, phasor in zip(color_vals, phasors):
            c = plt.cm.Set1(i)

            vl = Line2D([], [], color=c, linewidth=2)
            self.vec_lines.append(vl)
            ax_s.add_line(vl)

            vcl = Line2D([], [], color=c, linewidth=0.5)
            self.vec_connector_lines.append(vcl)
            ax_s.add_line(vcl)

            vd = Line2D([], [], color=c, marker='o', linewidth=0)
            self.vec_dots.append(vd)
            ax_s.add_line(vd)

            vtl = Line2D([], [], color=c)
            self.vec_trace_lines.append(vtl)
            ax_s.add_line(vtl)

            td = Line2D([], [], color=c, marker='o', linewidth=0)
            self.time_dots.append(td)
            ax_t.add_line(td)

            tcl = Line2D([], [], color=c, linewidth=0.5)
            self.time_connector_lines.append(tcl)
            ax_t.add_line(tcl)

            ttl = Line2D([], [], color=c)
            self.time_trace_lines.append(ttl)
            ax_t.add_line(ttl)

        animation.TimedAnimation.__init__(
            self, fig, interval=time_stretch/self.dt, blit=blit,
            repeat=repeat, repeat_delay=repeat_delay)

    def new_frame_seq(self):
        return iter(range(self.t.size))

    def _draw_frame(self, framedata):
        self._drawn_artists = []

        # advance phasors and plot them
        for i, phasor in enumerate(self.phasors):
            phasor.advance(self.dt)
            self.vec_lines[i].set_data([0, phasor.re], [0, phasor.im])
            self.vec_dots[i].set_data(phasor.re, phasor.im)
            self.vec_connector_lines[i].set_data([phasor.re, phasor.re],
                                                 [phasor.im, self.im_range[0]])
            self.vec_trace_lines[i].set_data(phasor.trace_re, phasor.trace_im)
            self.time_dots[i].set_data(phasor.re, phasor.t)
            self.time_connector_lines[i].set_data([phasor.re, phasor.re],
                                                  [0, phasor.t])
            self.time_trace_lines[i].set_data(phasor.trace_re, phasor.trace_t)

        # add lines to _drawn_artists
        for phasor in self.phasors:
            self._drawn_artists.extend(self.vec_lines)
            self._drawn_artists.extend(self.vec_dots)
            self._drawn_artists.extend(self.vec_connector_lines)
            self._drawn_artists.extend(self.vec_trace_lines)
            self._drawn_artists.extend(self.time_dots)
            self._drawn_artists.extend(self.time_connector_lines)
            self._drawn_artists.extend(self.time_trace_lines)

    def _init_draw(self):
        # clear the phasor trajectories
        for phasor in self.phasors:
            phasor.clear()

        # reset line data
        if getattr(self, '_drawn_artists', None) is not None:
            for a in self._drawn_artists:
                a.set_data([], [])
