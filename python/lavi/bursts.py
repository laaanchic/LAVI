"""Burst detection based on Bursts_detection_WTPL_v1_0_0.m.

The implementation keeps MATLAB-style one-based sample and frequency indices in
``burst_dists`` for direct comparison with the original toolbox. NumPy arrays
are indexed zero-based internally.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import maximum_filter
from scipy.signal import butter, sosfiltfilt

from .burst_helpers import get_burst_durations, get_freq_span, show_current_time, wtpl_burst_durations


class TFRFunction(Protocol):
    def __call__(self, data: NDArray[np.float64], fs: float, foi: NDArray[np.float64], width: float) -> NDArray[np.complex128]: ...


class WTPLFunction(Protocol):
    def __call__(self, fourier: NDArray[np.complex128], fs: float, foi: NDArray[np.float64], toi: NDArray[np.int64]) -> NDArray[np.float64]: ...


@dataclass
class BurstParameters:
    foi: NDArray[np.float64] | None = None
    wave_width: float = 5.0
    ptl_high: float = 90.0
    ptl_dur: float = 75.0
    SD_dur: float = 0.0
    downsmpFact: int = 1
    minQuietFs: float = 0.25
    output_power: bool = False
    art_thresh: float = np.nan
    lastSample: int | None = None  # MATLAB-style one-based sample
    min_f_gap: NDArray[np.float64] | None = None
    max_span: float | None = None
    wtpl_thr: float = 0.5
    ctrl_bursts: str = "no"
    fs: float | None = None
    min_dur: float = 1.0

    def resolved(self, n_times: int, fs: float) -> "BurstParameters":
        foi = np.asarray(self.foi if self.foi is not None else 10.0 ** np.arange(0.5, 1.6500001, 0.025), dtype=float)
        min_f_gap = np.asarray(self.min_f_gap if self.min_f_gap is not None else np.maximum(np.ceil(foi / 4.0), 4.0), dtype=float)
        return BurstParameters(
            **{
                **asdict(self),
                "foi": foi,
                "min_f_gap": min_f_gap,
                "max_span": float(self.max_span if self.max_span is not None else foi[-1] - foi[0]),
                "lastSample": int(self.lastSample if self.lastSample is not None else n_times),
                "fs": float(self.fs if self.fs is not None else fs),
            }
        )


@dataclass
class BurstOutput:
    pmtr: BurstParameters
    burst_dists: NDArray[np.float64]
    trough_peak: list[tuple[NDArray[np.int64], NDArray[np.int64]]]
    peak_span: list[NDArray[np.float64]]
    arts: NDArray[np.bool_]
    col_names: tuple[str, ...]
    pwr: NDArray[np.float64] | None = None
    filtd: NDArray[np.float64] | None = None


COL_NAMES = (
    "maxima freq (Hz)", "maxima timing (sample)", "begin (sample)", "end (sample)",
    "nearest trough", "duration (period)", "duration (ms)", "power (rel pctl)",
    "power (absolute)", "number of troughs", "number of peaks",
    "begin high percentile (sample)", "end high percentile (sample)", "ignore",
    "low Span freq (Hz)", "high Span freq (Hz)", "low Peak freq (Hz)",
    "high Peak freq (Hz)", "begin (sample, wtpl)", "end (sample, wtpl)",
    "duration WTPL (periods)",
)


def _regional_maxima(x: NDArray[np.float64]) -> NDArray[np.bool_]:
    """Approximate MATLAB ``imregionalmax`` for a 2-D matrix."""
    neighbourhood_max = maximum_filter(x, size=3, mode="constant", cval=-np.inf)
    return (x == neighbourhood_max) & (x > 0)


def _bandpass_rows(data: NDArray[np.float64], fs: float, foi: NDArray[np.float64]) -> NDArray[np.float64]:
    out = np.empty((foi.size, data.size), dtype=float)
    nyquist = fs / 2.0
    for k, freq in enumerate(foi):
        low, high = max(0.01, freq - 2.0), min(nyquist * 0.999, freq + 2.0)
        if low >= high:
            raise ValueError(f"Invalid band around {freq:g} Hz for fs={fs:g}")
        sos = butter(4, [low, high], btype="bandpass", fs=fs, output="sos")
        out[k] = sosfiltfilt(sos, data)
    return out


def detect_bursts_wtpl(
    data: ArrayLike,
    fs: float,
    *,
    pmtr: BurstParameters | None = None,
    tfr_func: TFRFunction,
    wtpl_func: WTPLFunction,
    artifact_func: Callable[[NDArray[np.float64], NDArray[np.float64], float, float], tuple[NDArray[np.float64], NDArray[np.bool_]]] | None = None,
    verbose: bool = True,
) -> tuple[BurstOutput, BurstParameters]:
    """Detect bursts and calculate burst statistics.

    ``tfr_func`` should return complex Fourier coefficients with shape
    ``(n_frequencies, n_times)``. ``wtpl_func`` should return WTPL with the
    same shape. These callbacks intentionally expose the two external MATLAB
    dependencies rather than replacing them with a non-equivalent algorithm.
    """
    raw = np.asarray(data, dtype=float).squeeze()
    if raw.ndim != 1:
        raise ValueError("data must contain one channel and one continuous trial")
    if not np.isfinite(fs) or fs <= 0:
        raise ValueError("fs must be positive and finite")

    p = (pmtr or BurstParameters()).resolved(raw.size, fs)
    foi = np.asarray(p.foi, dtype=float)
    n_f = foi.size
    if n_f < 2:
        raise ValueError("At least two frequencies are required")

    f1 = max(1.0, foi[0] - (foi[1] - foi[0]))
    f2 = foi[-1] + (foi[-1] - foi[-2])
    foi2 = np.concatenate(([f1], foi, [f2]))

    if verbose:
        print(f"Started calculating power at {show_current_time()}")
    fourier = np.asarray(tfr_func(raw[np.newaxis, :], p.fs, foi2, p.wave_width)).squeeze()
    if fourier.shape != (foi2.size, raw.size):
        raise ValueError(f"tfr_func returned {fourier.shape}; expected {(foi2.size, raw.size)}")
    offline_pwr = np.abs(fourier) ** 2

    margin = int(np.ceil(p.fs / f1)) + 1
    toi = np.arange(margin - 1, raw.size - margin, dtype=np.int64)
    WTPL = np.asarray(wtpl_func(fourier, p.fs, foi2, toi), dtype=float)
    if WTPL.shape != offline_pwr.shape:
        raise ValueError("wtpl_func must return the same frequency-by-time shape as tfr_func")

    if np.isfinite(p.art_thresh):
        if artifact_func is None:
            raise ValueError("artifact_func is required when art_thresh is finite")
        pwr_clean, arts = artifact_func(raw, offline_pwr, p.fs / p.downsmpFact, p.art_thresh)
    else:
        pwr_clean = offline_pwr.copy()
        arts = np.zeros(raw.size, dtype=bool)

    high = np.nanpercentile(pwr_clean, p.ptl_high, axis=1)
    dur = np.nanpercentile(pwr_clean, p.ptl_dur, axis=1) + p.SD_dur * np.nanstd(pwr_clean, axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        pwr_masked = pwr_clean / high[:, None]
        pwr_for_dur = pwr_clean / dur[:, None]
    pwr_masked[(pwr_masked < 1) | ~np.isfinite(pwr_masked)] = 0
    pwr_for_dur[(pwr_for_dur < 1) | ~np.isfinite(pwr_for_dur)] = 0

    bw = _regional_maxima(pwr_masked)
    fi_all, ti = np.nonzero(bw)
    keep = (fi_all != 0) & (fi_all != foi2.size - 1) & (ti < int(p.lastSample))
    fi = fi_all[keep] - 1
    ti = ti[keep]

    pwr_clean = pwr_clean[1:-1]
    pwr_masked = pwr_masked[1:-1]
    pwr_for_dur = pwr_for_dur[1:-1]
    WTPL = WTPL[1:-1]
    threshold = np.nanpercentile(WTPL, 75)
    wtpl = WTPL.copy()
    wtpl[wtpl < threshold] = 0

    if verbose:
        print(f"Started filtering at {show_current_time()}")
    filtd = _bandpass_rows(raw, p.fs, foi)

    n_bursts = ti.size
    burst = np.zeros((n_bursts, 21), dtype=float)
    trough_peak: list[tuple[NDArray[np.int64], NDArray[np.int64]]] = []
    peak_span: list[NDArray[np.float64]] = []
    too_short = np.zeros(n_bursts, dtype=bool)

    for x in range(n_bursts):
        burst[x, 0] = foi[fi[x]]
        burst[x, 13] = fi[x] + 1  # preserve MATLAB one-based frequency index
        burst[x, 1] = (ti[x] + 1) * p.downsmpFact

        _, _, _, _ = get_burst_durations(fi, ti, x, p, filtd, pwr_masked, filtered_freq_indices=fi)
        t_before, t_after, peaks, troughs = get_burst_durations(fi, ti, x, p, filtd, pwr_for_dur, filtered_freq_indices=fi)
        burst[x, 2:4] = ((t_before + 1) * p.downsmpFact, (t_after + 1) * p.downsmpFact)
        burst[x, 5] = (t_after - t_before) / np.ceil(p.fs / burst[x, 0])
        too_short[x] = burst[x, 5] < p.min_dur

        if not too_short[x]:
            burst[x, 6] = (t_after - t_before) / p.fs * 1000.0
            burst[x, 7] = pwr_masked[fi[x], ti[x]]
            burst[x, 8] = pwr_clean[fi[x], ti[x]]
            trough_peak.append((troughs + 1, peaks + 1))
            burst[x, 9:11] = (troughs.size, peaks.size)

            radius = int(np.ceil(p.fs / burst[x, 0] * 1.5))
            lo, hi = max(1, ti[x] - radius), min(raw.size - 2, ti[x] + radius)
            candidates = np.arange(lo, hi + 1)
            is_trough = (filtd[fi[x], candidates] < filtd[fi[x], candidates - 1]) & (filtd[fi[x], candidates] < filtd[fi[x], candidates + 1])
            locs = candidates[is_trough]
            if locs.size:
                distance = np.abs(locs - ti[x])
                nearest = locs[distance == distance.min()]
                if nearest.size > 1:
                    nearest = nearest[np.argmax(np.abs(filtd[fi[x], nearest])):np.argmax(np.abs(filtd[fi[x], nearest])) + 1]
                burst[x, 4] = nearest[0] + 1
            else:
                burst[x, 4] = ti[x] + 1
        else:
            trough_peak.append((np.empty(0, dtype=int), np.empty(0, dtype=int)))

        _, ind_above, ind_below = get_freq_span(pwr_masked[:, ti[x]], fi[x])
        burst[x, 14] = foi[fi[x] - ind_below]
        burst[x, 15] = foi[fi[x] + ind_above]
        f_ind = np.arange(fi[x] - ind_below, fi[x] + ind_above + 1)

        if f_ind.size > 1:
            start, stop = int(burst[x, 2]) - 1, int(burst[x, 3])
            curr = pwr_masked[f_ind, start:stop]
            maxima = np.max(curr, axis=0)
            argmax = np.argmax(curr, axis=0)
            argmax = argmax[maxima != 0]
            peak_f_ind = f_ind[argmax]
            if peak_f_ind.size:
                burst[x, 16] = foi[f_ind[argmax.min()]]
                burst[x, 17] = foi[f_ind[argmax.max()]]
                peak_span.append(foi[peak_f_ind])
            else:
                peak_f_ind = np.asarray([fi[x]])
                burst[x, 16:18] = foi[fi[x]]
                peak_span.append(foi[peak_f_ind])
        else:
            peak_f_ind = np.asarray([fi[x]])
            burst[x, 16:18] = foi[fi[x]]
            peak_span.append(foi[peak_f_ind])

        tb, ta = wtpl_burst_durations(peak_f_ind, burst[x], pwr_for_dur)
        burst[x, 11:13] = ((tb + 1) * p.downsmpFact, (ta + 1) * p.downsmpFact)
        burst[x, 5] = (ta - tb) / np.ceil(p.fs / burst[x, 0])
        tb, ta = wtpl_burst_durations(peak_f_ind, burst[x], wtpl)
        burst[x, 18:20] = (tb + 1, ta + 1)
        burst[x, 20] = (ta - tb) / np.ceil(p.fs / burst[x, 0])

    too_wide = (burst[:, 15] - burst[:, 14]) > float(p.max_span)
    keep = ~(too_short | too_wide)
    bst = burst[keep].copy()
    trough_peak = [v for v, k in zip(trough_peak, keep) if k]
    peak_span = [v for v, k in zip(peak_span, keep) if k]

    # Pruning. This intentionally preserves the MATLAB v1.0.0 energy formula,
    # including e2's (end + begin) expression, for reproducibility.
    with np.errstate(divide="ignore", invalid="ignore"):
        pwr = pwr_clean / high[1:-1, None]
    for bi in range(len(bst)):
        if bst[bi, 0] <= 0:
            continue
        overlap = ((bst[:, 2] >= bst[bi, 2]) & (bst[:, 2] <= bst[bi, 3])) | ((bst[:, 3] <= bst[bi, 3]) & (bst[:, 3] >= bst[bi, 2])) | ((bst[:, 2] <= bst[bi, 2]) & (bst[:, 3] >= bst[bi, 3]))
        overlap[bi] = False
        for ci in np.flatnonzero(overlap):
            idx1, idx2 = int(bst[bi, 13]) - 1, int(bst[ci, 13]) - 1
            if abs(bst[bi, 0] - bst[ci, 0]) <= np.asarray(p.min_f_gap)[idx1]:
                pwr1 = pwr[idx1, int(bst[bi, 1]) - 1]
                pwr2 = pwr[idx2, int(bst[ci, 1]) - 1]
                e1 = pwr1 * (bst[bi, 3] - bst[bi, 2])
                e2 = pwr2 * (bst[ci, 3] + bst[ci, 2])
                if e2 > e1:
                    bst[ci, 2] = min(bst[ci, 2], bst[bi, 2])
                    bst[ci, 3] = max(bst[ci, 3], bst[bi, 3])
                    bst[bi] = 0
                    break
                bst[bi, 2] = min(bst[ci, 2], bst[bi, 2])
                bst[bi, 3] = max(bst[ci, 3], bst[bi, 3])
                bst[ci] = 0

    final_keep = bst[:, 0] != 0
    bst = bst[final_keep]
    trough_peak = [v for v, k in zip(trough_peak, final_keep) if k]
    peak_span = [v for v, k in zip(peak_span, final_keep) if k]

    output = BurstOutput(
        pmtr=p,
        burst_dists=bst,
        trough_peak=trough_peak,
        peak_span=peak_span,
        arts=np.asarray(arts, dtype=bool),
        col_names=COL_NAMES,
        pwr=(pwr_clean / high[1:-1, None]) if p.output_power else None,
        filtd=filtd if p.output_power else None,
    )
    return output, p
