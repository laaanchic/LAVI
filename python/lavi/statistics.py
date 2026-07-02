"""Utilities for generating/significance limits.

``general_sig_lims`` is a Pythonic function version of the script
``General_SigLims.m``. It can be slow, just like the MATLAB script.
"""

from __future__ import annotations

import numpy as np
from scipy.io import savemat

from .core import prepare_lavi
from .spectra import get_pink_iafft_coefs_random_ap
from .surrogates import iaaft_loop_1d


def general_sig_lims(
    output_file: str | None = None,
    fs_values=(250, 500, 1000, 2000, 4000),
    dur_values=(120, 180, 300, 600),
    b_values=None,
    a: float = 1.0,
    reps: int = 20,
    foi=None,
    width: float = 5.0,
    lag: float = 1.5,
    max_iterate: int = 1000,
) -> tuple[np.ndarray, dict]:
    """Generate a lookup table of LAVI significance limits.

    Returns
    -------
    siglim : ndarray, shape (dur, fs, b, freq, 2)
    pmtrsig : dict
    """
    fs_values = np.asarray(fs_values, dtype=float)
    dur_values = np.asarray(dur_values, dtype=float)
    if b_values is None:
        b_values = np.arange(-1.2, 0.0 + 1e-12, 0.4)
    b_values = np.asarray(b_values, dtype=float)
    if foi is None:
        foi = 10 ** np.arange(0.5, 1.65 + 1e-12, 0.025)
    foi = np.asarray(foi, dtype=float)

    siglim = np.zeros((len(dur_values), len(fs_values), len(b_values), len(foi), 2))
    for di, dur in enumerate(dur_values):
        for fsi, fs in enumerate(fs_values):
            n = int(round(fs * dur))
            for bi, b in enumerate(b_values):
                coefs, _ = get_pink_iafft_coefs_random_ap(n, fs, foi, a, b)
                lavi_reps = np.zeros((reps, len(foi)))
                for ri in range(reps):
                    pink, _, _, _ = iaaft_loop_1d(coefs, np.sort(np.random.rand(n)), max_iterate)
                    cfg = {"foi": foi, "fs": fs, "lag": lag, "width": width, "verbose": False}
                    lavi_reps[ri, :] = prepare_lavi(cfg, pink)[0].squeeze()
                siglim[di, fsi, bi, :, :] = np.stack([np.min(lavi_reps, axis=0), np.max(lavi_reps, axis=0)], axis=1)

    pmtrsig = {
        "B": b_values,
        "DUR": dur_values,
        "FS": fs_values,
        "f": foi,
        "dimord": "dur_fs_b_freq_min/max",
        "lag": lag,
        "width": width,
        "script": "general_sig_lims.py",
    }
    if output_file is not None:
        savemat(output_file, {"SIGLIM": siglim, "pmtrSIG": pmtrsig})
    return siglim, pmtrsig
