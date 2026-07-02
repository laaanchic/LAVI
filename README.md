# LAVI
Toolbox to compute the rhythmicity profile and automatic band detection (both sustained and transient bands) of electrophysiological data, introduced in [Universal rhythmic architecture uncovers two modes of neural dynamics](https://www.nature.com/articles/s41467-026-73553-8) (Karvat et al 2026).

The Lagged-Angle Vector Index (LAVI) defines the rhythmicity spectrum based on phase persistence, that is, how well the phase (angle) in each frequency and timepoint can predict a future phase. Based on LAVI, the Automated Band-Border detection Algorithm (ABBA) automatically detects band peaks/troughs and borders, with statistical inference on the single-subject level.

## Content

This repository contains two implementations of LAVI. The original Matlab code is inside `matlab/`, whereas a Python-based version is inside `python/`. Both subfolders can be downloaded and installed independently. Please, find specific details of each implementation inside the corresponding subfolders.

## References
When using this toolbox please cite the following publications:
1. Karvat, G., Crespo-García, M., Vishne, G., Anderson, M. & Landau, A. N. (2026). Universal rhythmic architecture uncovers two modes of neural dynamics. Nature Communications. https://doi.org/10.1038/s41467-026-73553-8.
2. Oostenveld, R., Fries, P., Maris, E. & Schoffelen, J.-M. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Intell. Neuroscience 2011, 1:1-1:9 (2011).
3. Venema, V., Ament, F. & Simmer, C. A Stochastic Iterative Amplitude Adjusted Fourier Transform algorithm with improved accuracy. Nonlinear Processes in Geophysics 13, 321–328 (2006).
