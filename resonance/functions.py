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
