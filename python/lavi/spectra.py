"""Spectral helper functions for LAVI surrogate generation."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import minimize
from scipy.signal import welch, windows


def get_fft_frequencies(fs: float, n: int) -> np.ndarray:
    """Translate ``getFrequenciesOfFFT.m`` exactly.

    Unlike ``np.fft.fftfreq`` for even ``n``, this keeps the Nyquist bin as
    ``+fs/2`` and makes the following bins negative.
    """
    f = np.arange(n, dtype=float) * fs / n
    if n % 2:  # odd
        f[(n + 1) // 2 :] -= fs
    else:  # even; MATLAB changes N/2+2:end, i.e. 0-based N/2+1 onward
        f[n // 2 + 1 :] -= fs
    return f


def pwelch_to_amplitude(pxx: ArrayLike, f: ArrayLike, w: ArrayLike) -> np.ndarray:
    """Convert Welch PSD to amplitude, translated from ``Pwelch2amplitude.m``."""
    pxx = np.asarray(pxx, dtype=float)
    f = np.asarray(f, dtype=float)
    w = np.asarray(w, dtype=float)
    fbin = f[1] - f[0]
    cg = np.sum(w) / len(w)
    ng = np.sum(w**2) / len(w)
    return np.sqrt(pxx * (ng * fbin) / (cg**2) * 2.0)


def pwelch_nan(
    x: ArrayLike,
    window: ArrayLike,
    noverlap: int | None = None,
    nfft: int | None = None,
    fs: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Welch PSD while dropping any segment containing NaNs.

    Translation of ``pwelchNaN.m``. Input can be 1D or shape
    ``(n_times, n_channels)``.
    """
    X = np.asarray(x, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    window = np.asarray(window, dtype=float)
    L = len(window)
    if noverlap is None:
        noverlap = L // 2
    if nfft is None:
        nfft = int(2 ** np.ceil(np.log2(L)))
    step = int(L - noverlap)
    half_idx = int(np.floor(nfft / 2) + 1)
    PXX = np.zeros((half_idx, X.shape[1]))
    f = np.arange(nfft) * (fs / nfft)

    for chi in range(X.shape[1]):
        sig = X[:, chi]
        num_segments = int(np.floor((len(sig) - noverlap) / step))
        specs = []
        U = np.sum(window**2) / L
        for i in range(num_segments):
            start = i * step
            end = start + L
            if end > len(sig):
                break
            segment = sig[start:end] * window
            if np.isnan(segment).any():
                continue
            Xf = np.fft.fft(segment, nfft)
            specs.append((np.abs(Xf) ** 2) / (L * fs))
        if not specs:
            Pxx = np.full(nfft, np.nan)
        else:
            Pxx = np.sum(specs, axis=0) / (len(specs) * U)
        one = Pxx[:half_idx].copy()
        if len(one) > 2:
            one[1:-1] *= 2
        PXX[:, chi] = one
    return PXX.squeeze(), f[:half_idx]


def fit_gk(x: ArrayLike, y: ArrayLike) -> tuple[float, float]:
    """Fit ``y = a * x**b`` using Nelder-Mead, translated from ``fit_GK.m``."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x = x[mask]
    y = y[mask]

    def error(params: np.ndarray) -> float:
        a, b = params
        return float(np.sum((y - a * x**b) ** 2))

    res = minimize(error, np.array([1e5, -0.5], dtype=float), method="Nelder-Mead")
    return float(res.x[0]), float(res.x[1])


def get_ap_of_power(data: ArrayLike, fs: float, flim: tuple[float, float] = (5, 40)) -> tuple[float, float]:
    """Estimate aperiodic power offset/slope, translated from ``get_AP_of_Power.m``."""
    sig = np.asarray(data, dtype=float).squeeze()
    fs_i = int(round(fs))
    winsize = min(len(sig), 2 * fs_i)
    w = windows.hann(winsize, sym=False)
    pff, pxx = welch(sig, fs=fs_i, window=w, nperseg=winsize, noverlap=None, nfft=None)
    ind = (pff >= flim[0]) & (pff <= flim[-1])

    def obj(params: np.ndarray) -> float:
        a, b = params
        return float(np.linalg.norm(a * pff[ind] ** b - pxx[ind]))

    res = minimize(obj, np.array([1.0, 1.0]), method="Nelder-Mead")
    return float(res.x[0]), float(res.x[1])


def get_pink_iafft_coefs_pow(
    eeg: ArrayLike,
    window: ArrayLike,
    foi: ArrayLike,
    fs: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate IAAFT coefficients by fitting the observed aperiodic amplitude."""
    EEG = np.asarray(eeg, dtype=float).squeeze()
    w = np.asarray(window, dtype=float)
    foi = np.asarray(foi, dtype=float)
    N = len(EEG)
    if np.isnan(EEG).any():
        pxx, pff = pwelch_nan(EEG, w, None, None, fs)
    else:
        pff, pxx = welch(EEG, fs=fs, window=w, nperseg=len(w), noverlap=None, nfft=None)
    ff = get_fft_frequencies(fs, N)
    posf = ff[ff > 0]
    amp = pwelch_to_amplitude(pxx, pff, w)
    amp = amp * N / 2.0
    fit_fs = (pff >= foi[0]) & (pff <= foi[-1])
    a, b = fit_gk(pff[fit_fs], amp[fit_fs])
    ap = a * posf**b
    if N % 2:
        coefs = np.concatenate([[0.0], ap, ap[::-1]])
    else:
        coefs = np.concatenate([[0.0], ap, ap[:-1][::-1]])
    return coefs, amp


def get_pink_iafft_coefs_random_ap(
    n: int,
    fs: float,
    foi: ArrayLike,
    a: float,
    b: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate IAAFT coefficients from a specified random aperiodic power law.

    Translation of ``get_pink_iafft_coefs_random_ap.m``. MATLAB refits the
    derived amplitude spectrum with Curve Fitting Toolbox; here we use
    ``fit_gk`` for the same power-law form.
    """
    w = windows.hann(int(round(fs * 2)), sym=False)
    eeg = np.random.rand(int(n))
    pff, _ = welch(eeg, fs=fs, window=w, nperseg=len(w), noverlap=None, nfft=None)
    ff = get_fft_frequencies(fs, int(n))
    posf = ff[ff > 0]
    pow_ = a * pff**b
    amp = pwelch_to_amplitude(pow_, pff, w)
    amp = amp * int(n) / 2.0
    fit_fs = (pff >= np.asarray(foi)[0]) & (pff <= np.asarray(foi)[-1])
    aa, bb = fit_gk(pff[fit_fs], amp[fit_fs])
    ap = aa * posf**bb
    if int(n) % 2:
        coefs = np.concatenate([[0.0], ap, ap[::-1]])
    else:
        coefs = np.concatenate([[0.0], ap, ap[:-1][::-1]])
    return coefs, amp
