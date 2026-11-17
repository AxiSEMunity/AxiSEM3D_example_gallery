#!/usr/bin/env python3
"""Low-pass at the mesh period, as used for Fig. 8: zero-phase Butterworth of
order 4, applied forward and backward (scipy sosfiltfilt), corner at the mesh period (Mars 20 s,
Enceladus 2 s). Before filtering, the last TAPER seconds of each record are tapered to zero with a
half cosine, so that the end of the record (which cuts through ongoing motion) does not create a
filter transient; nothing else is changed."""
import numpy as np
from scipy.signal import butter, sosfiltfilt

MESH_PERIOD = {'mars': 20.0, 'enceladus': 2.0}
END_TAPER = {'mars': 100.0, 'enceladus': 20.0}      # seconds (five and ten mesh periods)


def lowpass(x, dt, period, taper_s):
    """x: (..., nt) records sampled at dt; returns the low-passed records (float64)."""
    x = np.asarray(x, dtype=np.float64)
    n = int(round(taper_s / dt))
    w = np.ones(x.shape[-1])
    if n > 1:
        w[-n:] = 0.5 * (1.0 + np.cos(np.linspace(0.0, np.pi, n)))
    sos = butter(4, 1.0 / period, btype='low', fs=1.0 / dt, output='sos')
    return sosfiltfilt(sos, x * w, axis=-1)
