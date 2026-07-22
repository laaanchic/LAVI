"""Python port of the MATLAB LAVI toolbox."""

from .core import compute_lavi, prepare_lavi, prepare_lavi_mne
from .abba import abba
from .surrogates import compute_pink_lavi
from .validation import compare_arrays, compare_complex_parts
from .spectra import (
    get_ap_of_power,
    get_fft_frequencies,
    get_pink_iafft_coefs_pow,
    get_pink_iafft_coefs_random_ap,
    pwelch_nan,
    pwelch_to_amplitude,
    fit_gk,
)

"""Burst detection additions for the Python LAVI toolbox."""

from .burst_helpers import (
    get_burst_durations,
    get_freq_span,
    show_current_time,
    wtpl_burst_durations,
)
from .bursts import BurstOutput, BurstParameters, detect_bursts_wtpl

__all__ = [
    "compute_lavi",
    "prepare_lavi",
    "prepare_lavi_mne",
    "abba",
    "compute_pink_lavi",
    "compare_arrays",
    "compare_complex_parts",
    "get_ap_of_power",
    "get_fft_frequencies",
    "get_pink_iafft_coefs_pow",
    "get_pink_iafft_coefs_random_ap",
    "pwelch_nan",
    "pwelch_to_amplitude",
    "fit_gk",
    "BurstOutput",
    "BurstParameters",
    "detect_bursts_wtpl",
    "get_burst_durations",
    "get_freq_span",
    "show_current_time",
    "wtpl_burst_durations",
]