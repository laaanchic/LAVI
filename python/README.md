# LAVI Python

Python-based version of the **LAVI toolbox**, which was programmed in Matlab. This package prioritises **numerical equivalence with the MATLAB implementation** over *pythonic* style. In several places the code deliberately mirrors MATLAB behaviour, including rounding, indexing conventions, padding, wavelet construction and edge handling. However, the code has a different organization as the original toolbox, but keeps the original names of the funtions.

The original **LAVI Toolbox** was created to compute the rhythmicity profile and automatic band detection (both sustained and transient bands) of electrophysiological data, introduced in [Universal rhythmic architecture uncovers distinct modes of neural dynamics](https://www.biorxiv.org/content/10.1101/2024.12.05.627113v1) (Karvat et al 2024).

This repository includes Jupyter notebooks with walkthrough examples that illustrate how to apply the functions to data matrices and MNE-Python data objects. They also contain validation steps on small example datasets, comparing the outputs between Matlab and Python implementations. In the validation tests, the band peaks and borders were identical, while the LAVI indexes differ from the 5th-6th decimal digit on. LAVI profiles on figures look the same.

## Installation

Clone or download the repository and move to the project directory:

```bash
cd lavi_python
```

Install the package:

```bash
pip install .
```

This installs the LAVI toolbox into your current Python environment together with its required dependencies.

Once installed, the package can be imported from Python:

```python
import lavi
```

### For developers

If you plan to modify the source code or contribute to the project, install the package in **editable mode** instead:

```bash
pip install -e .
```

This installs the package into your current Python environment while keeping it linked to the source code. Any changes made to the files in the `lavi/` directory are immediately available without reinstalling the package.

### Running the tests

To verify that the installation was successful, run the test suite from the project root:

```bash
pytest
```

or

```bash
python -m pytest
```

The tests are intended to verify that the package imports correctly and that the main functions execute without errors. Numerical validation against the original MATLAB implementation should be performed separately using the validation procedures described in the documentation.

## Basic usage with NumPy arrays

```python
from lavi import prepare_lavi, abba

# data shape: channels x time
lavi, cfg = prepare_lavi({"fs": 1000, "verbose": True}, data)
borders, var_names, sigvect = abba(lavi, cfg["foi"])
```

## Basic usage with MNE

```python
import mne
from lavi import prepare_lavi_mne

raw = mne.io.read_raw_fif("sample_raw.fif", preload=True)
lavi, cfg = prepare_lavi_mne(raw, picks="eeg")
```

MNE is optional and is only needed for `prepare_lavi_mne()` workflows.

## Modules

- `lavi.core`: `compute_lavi`, `prepare_lavi`, `prepare_lavi_mne`
- `lavi.wavelets`: `wavelet_light`, `tfr_light`, MATLAB-style helper functions
- `lavi.abba`: automatic band/border assignment
- `lavi.spectra`: Welch, FFT and aperiodic-spectrum helpers
- `lavi.surrogates`: pink-noise surrogate and IAAFT routines
- `lavi.statistics`: slow lookup-table generation equivalent to `General_SigLims.m`
- `lavi.io`: minimal MATLAB/FieldTrip-like loading helpers
- `lavi.plotting`: simple LAVI plotting helpers
- `lavi.validation`: optional array-comparison helpers for MATLAB/Python validation

## Notes and caveats

- The convolution step in `tfr_light()` intentionally uses `scipy.signal.convolve(..., mode="same", method="direct")`; it has not been replaced with an FFT or alternative alignment.
- Indices returned by `abba()` are 0-based by default. Use `matlab_indices=True` for MATLAB-style 1-based indices in the `BegI`, `EndI` and `PeakI` columns.
- `general_sig_lims()` can be very slow, as in the original MATLAB script.
