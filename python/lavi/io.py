"""Input/output helpers for LAVI."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.io import loadmat


def load_fieldtrip_like_mat(path: str, variable: str = "data") -> Any:
    """Load a MATLAB file and return a FieldTrip-like object/struct.

    This is deliberately minimal because MATLAB structs vary substantially.
    For MNE workflows, prefer reading your data with MNE and then calling
    ``prepare_lavi_mne``.
    """
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    return mat[variable]


def fieldtrip_trial_to_array(data: Any, picks: list[int] | np.ndarray | None = None) -> tuple[np.ndarray, float | None, list[str] | None]:
    """Extract ``trial``, ``fs`` and labels from a simple FieldTrip-like struct."""
    arr = np.asarray(data.trial)
    if picks is not None:
        arr = arr[picks, :]
    fs = getattr(data, "fs", None)
    labels = getattr(data, "label", None)
    if labels is not None and picks is not None:
        labels = [labels[i] for i in picks]
    return arr, fs, labels
