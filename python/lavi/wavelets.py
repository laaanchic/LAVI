"""Wavelet transforms used by LAVI.

This module contains Python translations of ``waveletLight.m`` and
``tfrLight.m``. The implementation intentionally prioritises numerical
agreement with the MATLAB LAVI toolbox over Pythonic optimisation.

Several helper functions reproduce MATLAB behaviours that differ subtly from
NumPy/SciPy defaults, such as rounding halves away from zero and using a
MATLAB-like colon operator for floating-point ranges.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from scipy.signal import convolve


def matlab_round(x: ArrayLike) -> np.ndarray:
    """Return MATLAB-style rounded values.

    MATLAB rounds half-integers away from zero, whereas NumPy uses bankers'
    rounding. Keeping this behaviour explicit is important for validation
    against the MATLAB toolbox.
    """
    x = np.asarray(x)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def matlab_colon(start: float, step: float, stop: float) -> np.ndarray:
    """Approximate MATLAB ``start:step:stop`` behaviour.

    The final value is included only when it is reached by the step sequence,
    with a small tolerance to reduce floating-point boundary artefacts.
    """
    n = int(np.floor((stop - start) / step + 1e-12))
    return start + step * np.arange(n + 1)


def nextpow2(x: float) -> int:
    """Return the exponent used by MATLAB ``nextpow2`` for positive scalars."""
    return int(np.ceil(np.log2(x)))


def tfr_light(
    data: ArrayLike,
    fsample: float,
    foi: ArrayLike | str,
    width: float | ArrayLike,
    verbose: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """MATLAB-faithful Python translation of ``tfrLight.m``.

    Parameters
    ----------
    data : array, shape (n_channels, n_times)
        Input data.
    fsample : float
        Sampling frequency in Hz.
    foi : array-like or "all"
        Frequencies of interest.
    width : float or array-like
        Wavelet width in cycles.
    verbose : bool
        Print progress messages.

    Returns
    -------
    spectrum : complex ndarray, shape (n_channels, n_freqs, n_times)
        Complex time-frequency spectrum.
    foi : ndarray
        Effective output frequencies after MATLAB-style frequency binning.
    toi : ndarray
        Output times in seconds.
    """

    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError("data must have shape (n_channels, n_times)")

    nchan, ndatsample = data.shape

    # MATLAB ``nanmean(data)`` operates along the first non-singleton
    # dimension. For channels x time data, this subtracts the channel mean at
    # each time point, not the temporal mean of each channel.
    data = data - np.nanmean(data, axis=0, keepdims=True)

    time = np.arange(1, ndatsample + 1, dtype=float) / fsample
    toi = time.copy()
    gwidth = 3

    dattime = ndatsample / fsample
    pad = 2 ** nextpow2(dattime)
    endnsample = int(matlab_round(pad * fsample))
    endtime = pad

    if endnsample < ndatsample:
        raise ValueError("The padding is shorter than the data.")

    if isinstance(foi, str):
        if foi != "all":
            raise ValueError('foi must be array-like or "all"')
        freqboilim = matlab_round(
            np.array([0, fsample / 2]) / (fsample / endnsample)
        ).astype(int) + 1
        freqboi = np.arange(freqboilim[0], freqboilim[1] + 1)
    else:
        foi_in = np.asarray(foi, dtype=float)
        freqboi = matlab_round(foi_in / (fsample / endnsample)).astype(int) + 1
        freqboi = np.unique(freqboi)

    foi_out = (freqboi - 1) / endtime

    if foi_out[0] == 0:
        foi_out = foi_out[1:]
        freqboi = freqboi[1:]

    nfreqoi = len(foi_out)
    offset = int(matlab_round(time[0] * fsample))
    toi = np.unique(matlab_round(toi * fsample) / fsample)

    # MATLAB timeboi is 1-based; convert only at the final array indexing step.
    timeboi = matlab_round(toi * fsample - offset).astype(int) + 1
    ntimeboi = len(timeboi)

    width = np.asarray(width, dtype=float)
    if width.ndim == 0 or width.size == 1:
        width = np.ones(nfreqoi) * float(width)
    else:
        width = width.ravel()
    if len(width) != nfreqoi:
        raise ValueError("width must be scalar or have one value per frequency")

    wavelets: list[np.ndarray] = []
    taplen = np.zeros(nfreqoi, dtype=int)

    for ifreq in range(nfreqoi):
        dt = 1 / fsample
        sf = foi_out[ifreq] / width[ifreq]
        st = 1 / (2 * np.pi * sf)

        toi2 = matlab_colon(-gwidth * st, dt, gwidth * st)
        A = 1 / np.sqrt(st * np.sqrt(np.pi))
        tap = A * np.exp(-(toi2**2) / (2 * st**2))

        acttapnumsmp = tap.size
        taplen[ifreq] = acttapnumsmp

        ind = matlab_colon(
            -(acttapnumsmp - 1) / 2,
            1,
            (acttapnumsmp - 1) / 2,
        ) * ((2 * np.pi / fsample) * foi_out[ifreq])

        wavelet = tap * np.cos(ind) + 1j * tap * np.sin(ind)
        wavelets.append(wavelet.astype(np.complex128))

    spectrum = np.full(
        (nchan, nfreqoi, ntimeboi),
        np.nan + 1j * np.nan,
        dtype=np.complex128,
    )

    for ifreq in range(nfreqoi):
        if verbose:
            print(f"frequency {ifreq + 1} ({foi_out[ifreq]:.2f} Hz)")

        nsamplefreqoi = taplen[ifreq]
        req_mask = (
            (timeboi >= (nsamplefreqoi / 2))
            & (timeboi < (ndatsample - (nsamplefreqoi / 2)))
        )
        reqtimeboiind = np.where(req_mask)[0]
        reqtimeboi = timeboi[reqtimeboiind]

        if reqtimeboi.size > 0:
            dum = np.zeros((nchan, reqtimeboi.size), dtype=np.complex128)
            for ichan in range(nchan):
                # Keep this direct convolution path unchanged: validation showed
                # small LAVI differences (~1e-5) and identical band borders.
                dumconv = convolve(
                    data[ichan, :],
                    wavelets[ifreq],
                    mode="same",
                    method="direct",
                )
                dum[ichan, :] = dumconv[reqtimeboi - 1]

            spectrum[:, ifreq, reqtimeboiind] = dum

    return spectrum, foi_out, toi


def wavelet_light(
    data: ArrayLike,
    fsample: float,
    foi: ArrayLike,
    width: float,
) -> np.ndarray:
    """MATLAB-faithful Python translation of ``waveletLight.m``.

    Parameters
    ----------
    data : array, shape (n_channels, n_times)
        Time-domain data, or already Fourier-transformed data if complex.
    fsample : float
        Sampling frequency in Hz.
    foi : array-like
        Frequencies of interest.
    width : float
        Wavelet width in cycles.

    Returns
    -------
    spectrum : complex ndarray, shape (n_channels, n_freqs, n_times)
        Complex time-frequency spectrum.
    """
    data = np.asarray(data)
    foi = np.asarray(foi, dtype=float)

    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.ndim != 2:
        raise ValueError("data must have shape (n_channels, n_times)")

    n_chan, n_time = data.shape
    n_freq = len(foi)
    gwidth = 3

    signal_freq = (
        np.fft.fft(data, axis=1)
        if np.isrealobj(data)
        else data.astype(np.complex128, copy=False)
    )

    wltspctrm: list[np.ndarray] = []
    taplen = np.zeros(n_freq, dtype=int)

    for fi in range(n_freq):
        dt = 1 / fsample
        sf = foi[fi] / width
        st = 1 / (2 * np.pi * sf)

        toi = matlab_colon(-gwidth * st, dt, gwidth * st)
        A = 1 / np.sqrt(st * np.sqrt(np.pi))
        tap = A * np.exp(-(toi**2) / (2 * st**2))

        acttapnumsmp = tap.size
        taplen[fi] = acttapnumsmp

        ins = int(np.ceil(n_time / 2) - np.floor(acttapnumsmp / 2))
        prezer = np.zeros(ins)
        pstzer_len = int(n_time - ((ins - 1) + acttapnumsmp) - 1)
        pstzer = np.zeros(pstzer_len)

        ind = matlab_colon(
            -(acttapnumsmp - 1) / 2,
            1,
            (acttapnumsmp - 1) / 2,
        ) * ((2 * np.pi / fsample) * foi[fi])

        wavelet = np.concatenate([prezer, tap * np.cos(ind), pstzer]) + 1j * np.concatenate(
            [prezer, tap * np.sin(ind), pstzer]
        )

        if wavelet.size != n_time:
            raise RuntimeError(
                f"Wavelet length mismatch at frequency {foi[fi]} Hz: "
                f"{wavelet.size} != {n_time}"
            )

        wltspctrm.append(np.fft.fft(wavelet).reshape(1, n_time))

    timeboi = np.arange(1, n_time + 1)
    spectrum = np.full(
        (n_chan, n_freq, n_time),
        np.nan + 1j * np.nan,
        dtype=np.complex128,
    )

    for fi in range(n_freq):
        nsamplefreqoi = taplen[fi]
        req_mask = (
            (timeboi >= (nsamplefreqoi / 2))
            & (timeboi < (n_time - (nsamplefreqoi / 2)))
        )
        reqtimeboiind = np.where(req_mask)[0]
        reqtimeboi = timeboi[reqtimeboiind]

        if reqtimeboi.size > 0:
            dum = np.fft.ifft(
                signal_freq * np.tile(wltspctrm[fi], (n_chan, 1)),
                axis=1,
            )
            dum = np.fft.fftshift(dum, axes=1)
            dum = dum * np.sqrt(2 / fsample)
            spectrum[:, fi, reqtimeboiind] = dum[:, reqtimeboi - 1]

    return spectrum
