"""Validation helpers for comparing MATLAB and Python LAVI outputs.

These helpers are intentionally independent of the public LAVI computation
pipeline. They are useful when checking MATLAB reference files during package
development, but no debug dataclass or validation mode is exposed publicly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


@dataclass(frozen=True)
class ArrayComparison:
    """Summary statistics for two arrays."""

    max_abs_diff: float
    mean_abs_diff: float
    max_rel_diff: float
    correlation: float


def compare_arrays(reference: ArrayLike, candidate: ArrayLike) -> ArrayComparison:
    """Compare two numeric arrays after flattening finite entries.

    Complex arrays are compared using their absolute complex difference. The
    correlation is computed on real-valued flattened arrays; for complex values,
    use `compare_complex_parts` to inspect real and imaginary components
    separately.
    """
    ref = np.asarray(reference)
    cand = np.asarray(candidate)
    if ref.shape != cand.shape:
        raise ValueError(f"shape mismatch: {ref.shape} != {cand.shape}")

    diff = cand - ref
    abs_diff = np.abs(diff)
    denom = np.maximum(np.abs(ref), np.finfo(float).eps)
    rel_diff = abs_diff / denom

    finite = np.isfinite(np.asarray(ref).real) & np.isfinite(np.asarray(cand).real)
    finite = finite.ravel()
    ref_flat = np.asarray(ref).ravel()[finite]
    cand_flat = np.asarray(cand).ravel()[finite]

    if ref_flat.size > 1:
        correlation = float(np.corrcoef(np.real(ref_flat), np.real(cand_flat))[0, 1])
    else:
        correlation = np.nan

    return ArrayComparison(
        max_abs_diff=float(np.nanmax(abs_diff)),
        mean_abs_diff=float(np.nanmean(abs_diff)),
        max_rel_diff=float(np.nanmax(rel_diff)),
        correlation=correlation,
    )


def compare_complex_parts(reference: ArrayLike, candidate: ArrayLike) -> dict[str, ArrayComparison]:
    """Compare real, imaginary and magnitude parts of complex arrays."""
    ref = np.asarray(reference)
    cand = np.asarray(candidate)
    return {
        "real": compare_arrays(np.real(ref), np.real(cand)),
        "imag": compare_arrays(np.imag(ref), np.imag(cand)),
        "abs": compare_arrays(np.abs(ref), np.abs(cand)),
    }
