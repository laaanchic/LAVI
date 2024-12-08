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
%% Plot LAVI
figure(646); clf; hold on
set(gcf,'position',[680 400 560 420]);
plot(foi, LAVI); set(gca,'xscale','log','xtick',[2:2:10,20:10:foi(end)]);
legend(data.label)
%% Calculate the borders using ABBA
% The function ABBA accepts up to 5 input arguments:
% 1. LAVI, as calculated above by Prepare_LAVI
% 2. foi (frequencies of interest), should be the same as used for calculating LAVI
% 3. alpha_range: the frequency range in which we expact to find alpha. The
%      band of the peak in this band will be assigned the index 0, bands
%      with lower frequency will be assigned with negative indices, and
%      above with positive indices. Default: [6 14]
% 4. sig_lim: controls the definition of significance levels per frequency.
%    There are three alternatives to define sig_lim:
%    a. Using a table based of pre-calculated pink-noise simulations 
%       (that can be generated by the script General_SigLims.m) that
%       matches the data's duration, sampling rate, and aperiodic slope.
%    b. By simulating pink-noise that matches the power spectrum of the
%       data.
%    c. leave empty. In this case, the bands will be found, but without
%       meaningful statistical inference.
% 5. perFreq: a boolean defining whether the siginificance level per
%    frequency (1), or as the minimum/ maximum over all frequencies (0,
%    default).
% =========================================================================
% Option A: use a previously-saved table (much faster compared to option B)
% SIGLIM has 5 dimensions. 
% 1. duration (session dependent) 
% 2. sampling rate (session dependent) 
% 3. aperiodic slope ("b", channel dependent)
% 4. frequency (should be the same as LAVI)
% 5. min/ max.
% The input to ABBA is an Nfreq x 2 matrix, that is, the last two
% domensions relevant to the session/ channel
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
%% Option B: Generate pink noise matching the data and calculate its LAVI
cfg  =[];
cfg.Pink_reps = Pink_reps; 
cfg.durs = durs; 
cfg.foi = foi;
PINK = computePinkLAVI(cfg,dat(choi,:)); % dimord: rep_freq_chan
pink = permute(PINK,[3,2,1]); % dimord: chan_freq_rep
sig_lim = cat(3,min(pink,[],3), max(pink,[],3)); % dimord: chan_freq_min/max
[borders2,~,sigVect2] = ABBA(LAVI, foi, alpha_range, sig_lim, 0);
%% Option C: no meaningful statistical inference 
[borders3,~,sigVect3] = ABBA(LAVI, foi, alpha_range);
%%
% Borders: 1 x Ncahn cell array, with a matrix of bands' details per channel 
%          (each rows is a band, the meaning of columns explained in col_names).
% col_names: 11x1 cell array with the names of variables in the borders matrix. In detail:
%            1. BegI: the Index of band Beginning
%            2. EndI: the Index of band Ending
%            3. PeakI: the Index of band Peak/trough (maximal/ minimal LAVI)
%            4. BegF: the Frequency of band Beginning
%            5. EndF: the Frequency of band Ending
%            6. PeakF: the Frequency of band Peak/trough (maximal/ minimal LAVI)
%            7. PeakLAVI: the LAVI value at the peak/ trough
%            8. PeakRel: the difference between PeakLAVI and the median
%            9. Dir: Direction of the PeakLAVI. 1= Peak ("sustained bands"), 
%               -1= Trough ("transient bands").
%            10.Rel_alpha: band identity. Alpha= 0. Lower frequency bands 
%               get negative integers, Higher frequency positive integers.
%            11.Sig: bool, Significantce relative to pink noise simulations
%               (1- the peak/ trough was higher/ lower than the pink noise
%               distribution, 0 otherwise).
% sigVect: 1 x Ncahn cell array, with 1 x Nfreq vector per channel of
%          statistical inference per frequency. 
%% Plot LAVI and ABBA
figure(456); clf; hold on
set(gcf,'position',[580 480 560 420]);
rows  =3;
cols = 3;
cachol = [83 174 244]/255; % blue colour
varod = [230 89 106]/255; % pink colour
yarok = [108 192 12]/255; % green colour

for chi = 1:length(choi)
    subplot(rows,cols,chi); hold on
    if ~isempty(PINK)
        plot(foi, squeeze(PINK(:,:,chi)),'color',ones(1,3)*0.5,'color',varod);
    end
    plot(foi, LAVI(chi,:),'k','linewidth',1.5)
    ylim([0.1 0.8])
    posind = sigVect2{chi}>0;
    negind = sigVect2{chi}<0;
    scatter(foi(posind), ones(1,sum(posind))*0.1,35,yarok,'fill','s');
    scatter(foi(negind), ones(1,sum(negind))*0.1,35,cachol,'fill','s');

    set(gca,'xscale','log','xtick',[2:2:10,20:10:foi(end)])
    title(num2str(chi))
end