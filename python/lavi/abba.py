"""ABBA: automatic band/border assignment for LAVI profiles."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

VAR_NAMES = [
    "BegI",
    "EndI",
    "PeakI",
    "BegF",
    "EndF",
    "PeakF",
    "PeakLAVI",
    "PeakRel",
    "Dir",
    "Rel_alpha",
    "Sig",
]


def _deal_with_zeros(x: np.ndarray) -> np.ndarray:
    out = x.copy()
    zs = np.where(out == 0)[0]
    for z in zs:
        if z == 0:
            out[z] = out[z + 1] / 10.0 if len(out) > 1 else np.nan
        else:
            out[z] = out[z - 1] / 10.0
    return out


def abba(
    lavi: ArrayLike,
    foi: ArrayLike,
    alpha_range: ArrayLike | None = None,
    siglim: ArrayLike | None = None,
    per_freq: bool = False,
    matlab_indices: bool = False,
) -> tuple[list[np.ndarray], list[str], list[np.ndarray]]:
    """Find LAVI bands, borders, and significance.

    Parameters mirror ``ABBA.m``. By default, returned indices are 0-based
    Python indices. Set ``matlab_indices=True`` to add 1 to the index columns.
    """
    LAVI = np.asarray(lavi, dtype=float)
    if LAVI.ndim == 1:
        LAVI = LAVI[np.newaxis, :]
    foi = np.asarray(foi, dtype=float)
    n_chan, n_freq = LAVI.shape

    if alpha_range is None:
        alpha = np.tile(np.array([6.0, 14.0]), (n_chan, 1))
    else:
        alpha = np.asarray(alpha_range, dtype=float)
        if alpha.ndim == 1:
            alpha = np.tile(alpha, (n_chan, 1))

    borders_all: list[np.ndarray] = []
    sigvect_all: list[np.ndarray] = []

    for ch in range(n_chan):
        lav = LAVI[ch, :]
        if siglim is None:
            sig_lim = np.tile(np.nanmedian(lav), (2, n_freq))
        else:
            sig_arr = np.asarray(siglim, dtype=float)
            if sig_arr.ndim == 3:
                sig_lim = np.squeeze(sig_arr[ch, :, :])
            elif sig_arr.size > 2:
                sig_lim = sig_arr
            else:
                sig_lim = np.tile(sig_arr.reshape(2, 1), (1, n_freq))
            if sig_lim.shape[0] != 2:
                sig_lim = sig_lim.T
            if sig_lim.shape[0] != 2:
                raise ValueError("Wrong definition of siglim")

        if siglim is None or not per_freq:
            sig_lim[0, :] = np.nanmin(sig_lim[0, :])
            sig_lim[1, :] = np.nanmax(sig_lim[1, :])

        sig_vect = np.zeros_like(foi, dtype=float)
        sig_vect[lav > sig_lim[1, :]] = 0.5
        sig_vect[lav < sig_lim[0, :]] = -0.5

        ref = np.nanmedian(lav)
        reref = _deal_with_zeros(lav - ref)
        siman = np.sign(reref)
        flipp = np.concatenate([np.diff(siman), [0]])
        flipp[np.isnan(flipp)] = 0
        border_points = np.where(flipp != 0)[0]
        starts = np.concatenate([[0], border_points + 1])
        ends = np.concatenate([border_points, [len(foi) - 1]])
        valid = starts <= len(foi) - 1
        starts = starts[valid]
        ends = ends[valid]

        borders = np.full((len(starts), 11), np.nan, dtype=float)
        borders[:, 0] = starts
        borders[:, 1] = ends
        borders[:, 8] = siman[ends]

        for bi in range(borders.shape[0]):
            inds = np.arange(int(borders[bi, 0]), int(borders[bi, 1]) + 1)
            peak_ind = inds[np.nanargmax(np.abs(reref[inds]))]
            borders[bi, 2] = peak_ind
            borders[bi, 6] = lav[peak_ind]
            borders[bi, 7] = reref[peak_ind]

        borders[:, 3:6] = np.round(foi[borders[:, 0:3].astype(int)], 1)
        w_inds = np.where(
            (borders[:, 5] <= alpha[ch, 1])
            & (borders[:, 5] >= alpha[ch, 0])
            & (borders[:, 8] > 0)
        )[0]
        if w_inds.size:
            a_ind = w_inds[np.nanargmax(borders[w_inds, 6])]
            borders[:, 9] = np.arange(borders.shape[0]) - a_ind
        else:
            borders[:, 3:6] = np.nan
            borders[:, 9] = np.nan

        peak_idx = borders[:, 2].astype(int)
        borders[:, 10] = (sig_vect[peak_idx] * borders[:, 8]) == 0.5
        for bi in np.where(~borders[:, 10].astype(bool))[0]:
            sig_vect[int(borders[bi, 0]) : int(borders[bi, 1]) + 1] = 0

        if matlab_indices:
            borders[:, 0:3] += 1
        borders_all.append(borders)
        sigvect_all.append(sig_vect)

    return borders_all, VAR_NAMES.copy(), sigvect_all
