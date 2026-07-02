"""Pink-noise surrogate generation for LAVI."""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from numpy.typing import ArrayLike
from scipy.signal import windows

from .core import compute_lavi, prepare_lavi
from .spectra import get_pink_iafft_coefs_pow
from .wavelets import wavelet_light


def iaaft_loop_1d(
    fourier_coeff: ArrayLike,
    sorted_values: ArrayLike,
    max_reps: int | None = 1000,
    error_threshold: float = 2e-4,
    speed_threshold: float = 1e-6,
) -> tuple[np.ndarray, int, float, float]:
    """IAAFT surrogate loop translated from ``iaaft_loop_1d.m``."""
    fourier_coeff = np.asarray(fourier_coeff)
    sorted_values = np.asarray(sorted_values, dtype=float).copy()
    if max_reps is None:
        max_reps = 1000
    error_amplitude = 1.0
    error_spec = 1.0
    old_total_error = 100.0
    speed = 1.0
    sd = np.std(sorted_values)
    if sd == 0 or not np.isfinite(sd):
        sd = 1.0

    idx = np.argsort(np.random.rand(*sorted_values.shape))
    y = np.empty_like(sorted_values, dtype=complex)
    y[idx] = sorted_values
    reps = 0

    while (
        (error_amplitude > error_threshold or error_spec > error_threshold)
        and speed > speed_threshold
        and reps < max_reps
    ):
        reps += 1
        old = y.copy()
        x = np.fft.ifft(y)
        phase = np.angle(x)
        x = fourier_coeff * np.exp(1j * phase)
        y = np.fft.fft(x)
        error_spec = float(np.mean(np.abs(np.real(y) - np.real(old))) / sd)

        old = y.copy()
        idx = np.argsort(np.real(y))
        y[idx] = sorted_values
        error_amplitude = float(np.mean(np.abs(np.real(y) - np.real(old))) / sd)

        total_error = error_spec + error_amplitude
        speed = abs((old_total_error - total_error) / total_error) if total_error != 0 else 0.0
        old_total_error = total_error

    return np.real(y), reps, error_amplitude, error_spec


def compute_pink_lavi(cfg: dict[str, Any] | None, data: ArrayLike) -> np.ndarray:
    """Generate pink-noise LAVI distributions, translated from ``computePinkLAVI.m``.

    Returns array with shape ``(pink_reps, n_freqs, n_channels)``.
    """
    cfg = {} if cfg is None else dict(cfg)
    pink_reps = int(cfg.get("Pink_reps", cfg.get("pink_reps", 20)))
    foi = np.asarray(cfg.get("foi", 10 ** np.arange(0.5, 1.65 + 1e-12, 0.025)), dtype=float)
    fs = float(cfg.get("fs", 1000))
    lag = float(cfg.get("lag", 1.5))
    width = float(cfg.get("width", 5))
    max_iterate = int(cfg.get("maxIterate", cfg.get("max_iterate", 1000)))

    arr = np.asarray(data, dtype=float)
    if arr.ndim == 1:
        arr = arr[np.newaxis, :]
    T = arr.shape[1] / fs
    durs = float(cfg.get("durs", T))

    if pink_reps == 0 or durs == 0:
        return np.empty((0, len(foi), arr.shape[0]))

    n_take = int(np.floor(min(durs, T) * fs))
    arr = arr[:, :n_take]
    n_chan, _ = arr.shape
    PINK = np.full((pink_reps, len(foi), n_chan), np.nan)
    w = windows.hann(int(round(fs * 2)), sym=False)
    t0 = time.time()

    for ch in range(n_chan):
        eeg = arr[ch, :]
        if np.isnan(eeg).any():
            eeg = eeg - np.nanmean(eeg)
        else:
            eeg = eeg - np.mean(eeg)
        coefs, _ = get_pink_iafft_coefs_pow(eeg, w, foi, fs)

        for ri in range(pink_reps):
            pink_noise, n_iter, _, _ = iaaft_loop_1d(coefs, np.sort(np.random.rand(*eeg.shape)), max_iterate)
            if cfg.get("verbose", True):
                print(
                    f"Running PINK ANALYSIS, Channel {ch + 1}/{n_chan} repeat {ri + 1}/{pink_reps} "
                    f"({n_iter} iterations). So far: {time.time() - t0:.2f} seconds."
                )
            for fi, freq in enumerate(foi):
                spectrum = wavelet_light(pink_noise, fs, np.array([freq]), width)[:, 0, :]
                PINK[ri, fi, ch] = compute_lavi(spectrum, fs, freq, lag)[0]
    return PINK
