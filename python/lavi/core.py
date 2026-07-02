"""Core LAVI functions.

The public functions in this module follow the MATLAB LAVI toolbox as closely
as possible. Numerical equivalence with MATLAB is prioritised over vectorised
or idiomatic alternatives when behaviours differ.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .wavelets import matlab_round, tfr_light, wavelet_light


@dataclass
class LAVIConfig:
    foi: np.ndarray | None = None
    fs: float = 1000.0
    lag: float = 1.5
    width: float = 5.0
    verbose: bool = True

    @classmethod
    def from_mapping(cls, cfg: dict[str, Any] | None = None, **kwargs: Any) -> "LAVIConfig":
        cfg = {} if cfg is None else dict(cfg)
        cfg.update(kwargs)
        foi = cfg.get("foi", None)
        if foi is None:
            foi = 10 ** np.arange(0.5, 1.65 + 1e-12, 0.025)
        return cls(
            foi=np.asarray(foi, dtype=float),
            fs=float(cfg.get("fs", cfg.get("sfreq", 1000.0))),
            lag=float(cfg.get("lag", 1.5)),
            width=float(cfg.get("width", 5.0)),
            verbose=bool(cfg.get("verbose", True)),
        )


def _as_2d(data: ArrayLike) -> np.ndarray:
    arr = np.asarray(data)
    if arr.ndim == 1:
        arr = arr[np.newaxis, :]
    if arr.ndim != 2:
        raise ValueError("data must have shape (n_channels, n_times) or (n_times,)")
    return arr

def compute_lavi(spectrum: ArrayLike, fs: float, f: float, lags: float = 1.5) -> np.ndarray:
    """Compute the Lagged Angle Vector Index for one frequency.

    Parameters
    ----------
    spectrum : complex array
        Shape ``(n_channels, n_times)`` or ``(n_channels, 1, n_times)``.
    fs : float
        Sampling frequency of the transformed data.
    f : float
        Frequency of interest in Hz.
    lags : float
        Lag in cycles. MATLAB default is 1.5.
    """
    spec = np.asarray(spectrum)
    if spec.ndim == 1:
        spec = spec[np.newaxis, :]
    elif spec.ndim == 3 and spec.shape[1] == 1:
        spec = spec[:, 0, :]
    elif spec.ndim != 2:
        raise ValueError("spectrum must have shape (n_channels, n_times) or (n_channels, 1, n_times)")

    n_ch, n_t = spec.shape
    # MATLAB uses round() with halves away from zero. This can affect
    # the lag sample width at .5 boundaries, so keep MATLAB-style rounding.
    width = int(matlab_round(lags / f * fs))
    if width <= 0 or width >= n_t:
        return np.full(n_ch, np.nan)

    sig0 = spec[:, : n_t - width]
    sig1 = spec[:, width:]

    # MATLAB removes NaNs based on channel 1 only.
    nanind = np.isnan(sig0[0, :]) | np.isnan(sig1[0, :])
    sig0 = sig0[:, ~nanind]
    sig1 = sig1[:, ~nanind]
    if sig0.shape[1] == 0:
        return np.full(n_ch, np.nan)

    a0 = np.abs(sig0)
    a1 = np.abs(sig1)
    numerator = np.sum(sig0 * np.conj(sig1), axis=1)
    denominator = np.sqrt(np.sum(a0 * a0, axis=1) * np.sum(a1 * a1, axis=1))
    out = np.abs(numerator / denominator)
    out[denominator == 0] = np.nan

    return out

def as_cell_array(items):
    out = np.empty((len(items),), dtype=object)
    for i, item in enumerate(items):
        out[i] = item
    return out

def prepare_lavi(cfg: dict[str, Any] | LAVIConfig | None, data: ArrayLike) -> tuple[np.ndarray, dict[str, Any]]:
    """Prepare data and compute LAVI across frequencies.

    This mirrors ``Prepare_LAVI.m`` but processes all channels together per
    frequency where possible.  Input data are ``channels x time``.
    """
    conf = cfg if isinstance(cfg, LAVIConfig) else LAVIConfig.from_mapping(cfg)
    arr = _as_2d(data)
    n_chan, _ = arr.shape
    foi = np.asarray(conf.foi, dtype=float)
    lavi = np.zeros((n_chan, len(foi)), dtype=float)
    t0 = time.time()

    for fi, freq in enumerate(foi):
        if conf.verbose:
            print(f"Running LAVI frequency {fi + 1}/{len(foi)} ({freq:.3g} Hz)")
        if np.isnan(arr).any():
            spectrum, foi_out, _ = tfr_light(arr, conf.fs, np.array([freq]), conf.width, verbose=False)
            # tfr_light may quantize the frequency based on padding.
            effective_freq = float(foi_out[0])
            spec = spectrum[:, 0, :]
        else:
            spec = wavelet_light(arr, conf.fs, np.array([freq]), conf.width)[:, 0, :]
            effective_freq = float(freq)
        
        lavi[:, fi] = compute_lavi(spec, conf.fs, effective_freq, conf.lag)
        
    out_cfg = {
        "foi": foi,
        "fs": conf.fs,
        "lag": conf.lag,
        "width": conf.width,
        "verbose": conf.verbose,
    }
    if conf.verbose:
        print(f"prepare_lavi took {time.time() - t0:.3g} seconds")
    return lavi, out_cfg


def prepare_lavi_mne(raw_or_epochs: Any, picks: Any = None, cfg: dict[str, Any] | None = None) -> tuple[np.ndarray, dict[str, Any]]:
    """Convenience wrapper for MNE Raw/Epochs objects.

    For Raw objects, this uses ``raw.get_data(picks=picks)``. For Epochs,
    epochs are concatenated along time. This is a minimal wrapper; project-
    specific preprocessing should be done before calling it.
    """
    if not hasattr(raw_or_epochs, "info") or not hasattr(raw_or_epochs, "get_data"):
        raise TypeError("Expected an MNE Raw/Epochs-like object with .info and .get_data().")
    data = raw_or_epochs.get_data(picks=picks)
    if data.ndim == 3:  # epochs x channels x times -> channels x concatenated time
        data = np.transpose(data, (1, 0, 2)).reshape(data.shape[1], -1)
    config = {} if cfg is None else dict(cfg)
    config.setdefault("fs", float(raw_or_epochs.info["sfreq"]))
    return prepare_lavi(config, data)
