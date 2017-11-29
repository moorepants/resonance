import numpy as np


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
