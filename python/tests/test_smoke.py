import numpy as np

from lavi import abba, compute_lavi, prepare_lavi
from lavi.wavelets import matlab_round, wavelet_light


def test_matlab_round_halves_away_from_zero():
    x = np.array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5])
    np.testing.assert_array_equal(matlab_round(x), np.array([-3, -2, -1, 1, 2, 3]))


def test_wavelet_and_lavi_smoke():
    rng = np.random.default_rng(0)
    data = rng.standard_normal((2, 1000))
    spec = wavelet_light(data, 1000, np.array([10.0]), 5.0)
    lavi = compute_lavi(spec[:, 0, :], 1000, 10.0)
    assert spec.shape == (2, 1, 1000)
    assert lavi.shape == (2,)


def test_prepare_lavi_and_abba_smoke():
    rng = np.random.default_rng(1)
    data = rng.standard_normal((2, 1000))
    foi = np.array([8.0, 10.0, 12.0])
    lavi, cfg = prepare_lavi({"fs": 1000, "foi": foi, "verbose": False}, data)
    borders, names, sigvect = abba(lavi, cfg["foi"])
    assert lavi.shape == (2, 3)
    assert len(borders) == 2
    assert names[0] == "BegI"
    assert len(sigvect) == 2
