clearvars
LAVIpath = 'LAVI_folder'; % give here the location in which the LAVI toolbox is saved
addpath(LAVIpath);
% data should contain a matrix with N_chan x N_timepoints preprocessed
% data, and the sampling frequency in Hz.
% The script can be tested with the provided 5 minutes rest session of 2 channel EEG
dataFileName = fullfile(LAVIpath, 'data'); % the path to saved data
siglimFileName = fullfile(LAVIpath, 'SIGLIM'); % the path to saved table of significance levels
load(dataFileName);
load(siglimFileName);
%
foi         = 10.^[0.5:0.025:1.65]; % frequencies of interest
fsample     = data.fs;              % sampling frequency
lag         = 1.5;                  % lag between the signal and its copy (in cycles, default = 1.5)
width       = 5;                    % wavelet width (in cycles, default = 5)
Pink_reps   = 20;                   % number of simulations created per channel. Default = 20.
durs        = 60;                   % the duration (in seconds) of each simulation. Default: the duration of the original signal
choi        = 1:3;                  % use this to limit the number of channels
%% Calculate LAVI of the data
cfg = [];
cfg.foi = foi;
cfg.fs = fsample;
cfg.lag = lag;
cfg.width = width;
cfg.verbose = 1;

dat = data.trial(choi,:);
if any(isnan(dat(:))) warning('Data contains NaNs. Calculating TFR by convolution in the time domain'); end 
[LAVI,cfg] = Prepare_LAVI(cfg,dat); % data has to be a matrix, channel x timepoints

%% Calculate the borders using ABBA
alpha_range = [6 14]; % the peak in LAVI in this range will be conidered alpha
siglim = zeros(size(dat,1), size(SIGLIM,4), 2);
dur = data.time(end)-data.time(1); % session duration is sec
[~,ind1] = min(abs(pmtrSIG.DUR - dur)); % the index of session duration in SIGLIM
[~,ind2] = min(abs(pmtrSIG.FS - data.fs)); % the index sampling frequency in SIGLIM
for ch = 1%:size(dat,1)
    [~,b] = get_AP_of_Power(dat(ch,:),data.fs,foi);
    [~,ind3] = min(abs(pmtrSIG.B - b)); % the index aperiodic slope in SIGLIM
    siglim(ch,:,:) = squeeze(SIGLIM(ind1,ind2,ind3,:,:));
end
[borders,col_names,sigVect] = ABBA(LAVI, foi, alpha_range, siglim, 0);

