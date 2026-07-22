"""Helper functions for burst detection translated from MATLAB."""

from __future__ import annotations

from datetime import datetime
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray


def show_current_time() -> str:
    """Return local time as ``H:MM:SS`` (translation of showCurrnetTime.m)."""
    now = datetime.now()
    return f"{now.hour}:{now.minute:02d}:{now.second:02d}"


def get_freq_span(col: ArrayLike, peak_index: int) -> tuple[int, int, int]:
    """Find contiguous positive frequency bins around a peak.

    Parameters
    ----------
    col
        One time-column of the thresholded power matrix.
    peak_index
        Zero-based frequency index of the burst peak.

    Returns
    -------
    freq_span, ind_above, ind_below
        Number of additional positive bins in total, above the peak, and below
        the peak. Names follow the MATLAB function, where ``above`` means
        increasing array index.
    """
    values = np.asarray(col)
    if values.ndim != 1:
        raise ValueError("col must be one-dimensional")
    if not 0 <= peak_index < values.size:
        raise IndexError("peak_index is outside col")

    ind_below = 0
    while peak_index - (ind_below + 1) >= 0:
        if values[peak_index - (ind_below + 1)] > 0:
            ind_below += 1
        else:
            break

    ind_above = 0
    while peak_index + (ind_above + 1) < values.size:
        if values[peak_index + (ind_above + 1)] > 0:
            ind_above += 1
        else:
            break

    return ind_above + ind_below, ind_above, ind_below


def _param(pmtr: Any, name: str) -> Any:
    if isinstance(pmtr, Mapping):
        return pmtr[name]
    return getattr(pmtr, name)


def get_burst_durations(
    freq_indices: ArrayLike,
    time_indices: ArrayLike,
    burst_index: int,
    pmtr: Any,
    filtered: ArrayLike,
    power_mask: ArrayLike,
    *,
    filtered_freq_indices: ArrayLike | None = None,
) -> tuple[int, int, NDArray[np.int64], NDArray[np.int64]]:
    """Determine burst boundaries and local peaks/troughs.

    This is the zero-based Python equivalent of ``get_burst_durations.m``.
    Returned sample indices are zero-based. The main wrapper converts them to
    MATLAB-style one-based samples in its output table.
    """
    fi = np.asarray(freq_indices, dtype=np.int64)
    ti = np.asarray(time_indices, dtype=np.int64)
    filtd = np.asarray(filtered)
    pwrr = np.asarray(power_mask)

    if filtered_freq_indices is None:
        fii = fi
    else:
        fii = np.asarray(filtered_freq_indices, dtype=np.int64)

    x = int(burst_index)
    f_idx = int(fi[x])
    filt_idx = int(fii[x])
    centre = int(ti[x])
    downsample_factor = int(_param(pmtr, "downsmpFact"))
    n = pwrr.shape[1]

    if centre <= 0 or centre >= filtd.shape[1] - 1:
        raise IndexError("burst peak must have neighbouring samples")

    peaks: list[int] = []
    troughs: list[int] = []

    if filtd[filt_idx, centre] > filtd[filt_idx, centre + 1] and filtd[filt_idx, centre] > filtd[filt_idx, centre - 1]:
        peaks.append(centre)
    if filtd[filt_idx, centre] < filtd[filt_idx, centre + 1] and filtd[filt_idx, centre] < filtd[filt_idx, centre - 1]:
        troughs.append(centre)

    t_before = centre * downsample_factor - 1
    while True:
        if t_before <= 0:
            break
        if filtd[filt_idx, t_before] > filtd[filt_idx, t_before + 1] and filtd[filt_idx, t_before] > filtd[filt_idx, t_before - 1]:
            peaks.append(t_before)
        if filtd[filt_idx, t_before] < filtd[filt_idx, t_before + 1] and filtd[filt_idx, t_before] < filtd[filt_idx, t_before - 1]:
            troughs.append(t_before)
        if not bool(pwrr[f_idx, t_before]):
            break
        t_before -= 1
    t_before += 1

    t_after = centre * downsample_factor + 1
    max_valid = min(n - 1, filtd.shape[1] - 2)
    while True:
        if t_after >= max_valid:
            break
        if filtd[filt_idx, t_after] > filtd[filt_idx, t_after + 1] and filtd[filt_idx, t_after] > filtd[filt_idx, t_after - 1]:
            peaks.append(t_after)
        if filtd[filt_idx, t_after] < filtd[filt_idx, t_after + 1] and filtd[filt_idx, t_after] < filtd[filt_idx, t_after - 1]:
            troughs.append(t_after)
        if not bool(pwrr[f_idx, t_after]):
            break
        t_after += 1
    t_after -= 1

    return (
        t_before,
        t_after,
        np.asarray(sorted(peaks), dtype=np.int64),
        np.asarray(sorted(troughs), dtype=np.int64),
    )


def wtpl_burst_durations(
    freq_indices: ArrayLike,
    burst_row: ArrayLike,
    wtpl: ArrayLike,
) -> tuple[int, int]:
    """Find burst boundaries across selected frequency bins.

    Parameters use zero-based indices internally. ``burst_row[1]`` is expected
    to contain a MATLAB-style one-based peak sample, as in ``burstDists``.
    Returned indices are zero-based.
    """
    f_ind = np.atleast_1d(np.asarray(freq_indices, dtype=np.int64))
    row = np.asarray(burst_row)
    arr = np.asarray(wtpl)
    if arr.ndim != 2:
        raise ValueError("wtpl must have shape (n_frequencies, n_times)")

    x0 = int(row[1]) - 1
    if not 0 <= x0 < arr.shape[1]:
        raise IndexError("burst peak sample is outside wtpl")

    t_before = x0 - 1
    while t_before > 0 and np.any(arr[f_ind, t_before]):
        t_before -= 1
    t_before += 1

    t_after = x0 + 1
    last = arr.shape[1] - 1
    while t_after < last and np.any(arr[f_ind, t_after]):
        t_after += 1
    t_after -= 1

    return t_before, t_after
