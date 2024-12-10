# timing
import numpy as np


def coherence(cs, ps1, ps2):
    """
    all arrays, cross spec, psds
    """
    coh = np.abs(cs)^2/(ps1*ps2)
    return coh

def phase_lag(cs_unnorm):
    return np.angle(cs_unnorm)

def time_lag(cs_unnorm):
    time_lag = phase_lag(cs_unnorm) / (2 * np.pi * cs_unnorm.freq)
    return time_lag
