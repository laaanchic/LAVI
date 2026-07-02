"""Small plotting helpers for LAVI outputs."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def plot_lavi(foi, lavi, labels=None, ax=None):
    """Plot LAVI profiles on a log-frequency x-axis."""
    foi = np.asarray(foi)
    LAVI = np.asarray(lavi)
    if LAVI.ndim == 1:
        LAVI = LAVI[np.newaxis, :]
    if ax is None:
        _, ax = plt.subplots()
    for ch in range(LAVI.shape[0]):
        label = None if labels is None else labels[ch]
        ax.plot(foi, LAVI[ch], label=label)
    ax.set_xscale("log")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("LAVI")
    if labels is not None:
        ax.legend()
    return ax
