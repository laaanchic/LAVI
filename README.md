# LAVI
Toolbox to compute the rhythmicity profile and automatic band detection (both sustained and transient bands) of electrophysiological data, introduced in [Universal rhythmic architecture uncovers distinct modes of neural dynamics](https://www.biorxiv.org/content/10.1101/2024.12.05.627113v1) (Karvat et al 2024).

The Lagged-Angle Vector Index (LAVI) defines the rhythmicity spectrum based on phase persistence, that is, how well the phase (angle) in each frequency and timepoint can predict a future phase. Based on LAVI, the Automated Band-Border detection Algorithm (ABBA) automatically detects band peaks/troughs and borders, with statistical inference on the single-subject level.

## Dependencies
Currently, LAVI and ABBA work on Matlab, with the Curve Fitting Toolbox.  
Although LAVI does not depend on the freely-available toolbox [Fieldtrip](https://www.fieldtriptoolbox.org/), some functions are inspired by Fieldtrip. Therefore, when using LAVI, please cite also the Fieldtrip toolbox (see reference below).

## Input
The main function, Prepare_LAVI, takes an N_channel x N_timepoints raw electrophysiological data (EEG, MEG, LFP, etc. Examplary EEG data available here as data.mat), and a configuration structure, including the frequencies of interest, sampling rate, lag duration and wavelt width.  
The data should be preprocessed according to conventional standards, that is, removed eye-blinks and 50/60 hz hum. Other artifacts, such as muscle artifacts, can be replaced with NaNs. We tested the algorhythm successfully with resting-state data with eyes closed without preprocessing. 

## How to use
Download the functions, and add the folder containing them to your Matlab path.  
The main functions to compute LAVI and ABBA are Prepare_LAVI.m and ABBA.m . Both functions contain help sections explaining their inuputs and outputs.  
For quick implenetation of the main functionalities of LAVI and ABBA, use LAVI_main.m .  
For an annotated, step-by-step, and detailed description of how to use LAVI and ABBA, see the script LAVI_walk_through.m (or live script LAVI_walk_through.mlx)

## References
When using this toolbox please cite the following publications:
1. Karvat, G., Crespo-García, M., Vishne, G., Anderson, M. & Landau, A. N. Universal rhythmic architecture uncovers distinct modes of neural dynamics. 2024.12.05.627113 Preprint at https://doi.org/10.1101/2024.12.05.627113 (2024).
2. Oostenveld, R., Fries, P., Maris, E. & Schoffelen, J.-M. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Intell. Neuroscience 2011, 1:1-1:9 (2011).
3. Venema, V., Ament, F. & Simmer, C. A Stochastic Iterative Amplitude Adjusted Fourier Transform algorithm with improved accuracy. Nonlinear Processes in Geophysics 13, 321–328 (2006).
