function [Output,pmtr] = Bursts_detection_WTPL_v1_0_0 (raw,pmtr)
% [Output,pmtr] = Bursts_detection_WTPL_v1_0_0 (raw,pmtr)
% Detects bursts and their statistics, as used in: 
% Karvat, G., Crespo-García, M., Vishne, G., Anderson, M. C., and 
% Landau, A. N. (2026). Universal rhythmic architecture uncovers two modes 
% of neural dynamics. Nature Communications. https://doi.org/10.1038/s41467-026-73553-8
%
% Dependencies: 
% ------------
% Fieldtrip toolbox: Available at https://www.fieldtriptoolbox.org/. 
%                    Please cite the original Fieldtrip paper when using 
%                    this function (doi: 10.1155/2011/156869)
% WTPL toolbox:      Available at https://zenodo.org/records/10664624. 
%                    If using WTPL for reporting burst duration, please cite the 
%                    original publication: https://doi.org/10.1162/jocn_a_02088
%
% Input: 
% -----
% raw:  a raw Fieldtrip structure, with one channel and one trial (of the whole session).
% pmtr: a structure with parameters. If empty, assigns default values.
% .foi:           frequencies of interest, in Hz. 
%                 Default: 10.^(0.5:0.025:1.65); 
% .wave_width:    width of the wavelet, in periods. Default: 5.
% .ptl_high:      threshold percentile for bursts. Default: 90.
% .ptl_dur:       percentile to use for determining duration. Default:75.
% .SD_dur:        value (in Standard Deviation) to add to the median for 
%                 duration calculation. Little et al 2019 used 1.75. Default:0.
% .downsmpFact:   downsampling parameter of the tfr. Default: 1 (no downsampling).
% .minQuietFs:    minimal of frequencies in foi below the thresold to be 
%                 considered "quiet". Default: 0.25.
% .output_power:  boolean, Keep the normalized power in output? Default: 0.
% .art_thresh:    optional- threshold for artifact removal in uV 
%                 (replacing 0.5 sec around artifacts with NaN). Skipped 
%                 if NaN. Default: NaN.
% .lastSample:    last sample to take into account. Default: end of session.
% .min_f_gap:     the minimum gap betwween overlapping bursts, in Hz. 
%                 Default: a quarter of each frequency, with a minimum of 4. 
% .max_span:      the maximum allowed frequency span. Default: the whole spectrum
% .wtpl_thr:      Threshold of WTPL to calculate duration. Default: 0.5.
% .ctrl_bursts:   'yes' or 'no', whether or not to plot individual bursts 
%                 for control. Default: 'no'.
% .fs:            sampling frequencies (in Hz). If not provided, tries 
%                 raw.fsample, and if not available, computes from raw.time{1};             % if fs (sapmling rate) not given, compute it  
% .min_dur        the minimum duration of a burst, in periods. Use 0.5 if 
%                 the duration is calculated according to the high percentile.
%                 Default: 1;                  % 
%
% Output: 
% ------
% .pmtr: structure with the updated parameter values (including defaults).
% .burstDists: N_burst x 21 matrix with columns:
%       1: maxima freq(Hz), 2: maxima timing(sample), 3: begin (sample), 
%       4: end (sample), 5: nearest trough, 6: duration (periods), 
%       7: duration (ms), 8: power (rel pctl), 9: power (absolute), 
%       10: number of troughs, 11: number of peaks, 
%       12: begin high percentile (sample), 13: end high percentile (sample)
%       14: freq index (ignore), 15: lowest SPAN frequency (Hz), 
%       16: highest SPAN frequency (Hz), 17: lower PEAK freq (Hz), 
%       18: highest PEAK freq (Hz), 19: begin (sample, wtpl), 
%       20: end(sample, wtpl), 21: duration WTPL (periods)
% .TroughPeak:  N_burst x 2 cell array, with the samples of all troughs (:,1) 
%               and all peaks (:.2) within the burst.
% .PeakSpan:    N_burst x 1 cell array, with the dominant (highest power)
%               frequency in each sample of the burst.
% .arts:        N_time x 1 (logical). Samples with detected artifacts 
%               (relevant when pmtr.arts is provided)
% .col_names:   21 x 1 cell array, names of the columns of the burstDists
%               matrix.
% .pwr:         N_foi x N_time double, the normalized time-frequency-represntation
%               (given if pmtr.output_power is true).
% .filtd:       N_foi x N_time double, the bandpass-filtered data (given if 
%               pmtr.output_power is true).
%
% Version 1.0.0, 19/06/2026 (Golan Karvat)
%
% To do: define the duration of the burst based on lagged-coherence,
% as suggested by Schmidt et al 2023, https://www.sciencedirect.com/science/article/pii/S0959438823001216

t_func = tic;
% get the parameters
if nargin<2; pmtr=[]; end

% Define parameters
if ~isfield(pmtr,'foi');         pmtr.foi = 10.^(0.5:0.025:1.65); end  % frequencies of interest, in Hz
if ~isfield(pmtr,'wave_width');  pmtr.wave_width   = 5; end            % the width of the wavelet, in periods
if ~isfield(pmtr,'ptl_high');    pmtr.ptl_high     = 90; end           % the threshold percentile for bursts
if ~isfield(pmtr,'ptl_dur');     pmtr.ptl_dur      = 75; end           % the percentile to use for determining duration
if ~isfield(pmtr,'SD_dur');      pmtr.SD_dur       = 0; end            % value (in SD) to add to the median for duration calculation. Little et al 2019 used 1.75
if ~isfield(pmtr,'downsmpFact'); pmtr.downsmpFact  = 1; end            % the downsampling parameter of the tfr
if ~isfield(pmtr,'minQuietFs');  pmtr.minQuietFs   = 0.25; end         % minimal of frequencies in foi below the thresold to be considered "quiet"
if ~isfield(pmtr,'output_power');pmtr.output_power = 0; end            % Keep the normalized power in output? 
if ~isfield(pmtr,'art_thresh');  pmtr.art_thresh   = NaN; end          % optional- threshold for artifact removal in uV (replacing 0.5 sec around artifacts with NaN). Skipped if NaN
if ~isfield(pmtr,'lastSample');  pmtr.lastSample   = raw.sampleinfo(end); end %last sample to take into account. Default: end of session
if ~isfield(pmtr,'min_f_gap');   pmtr.min_f_gap    = ceil(pmtr.foi./4); pmtr.min_f_gap(pmtr.min_f_gap<4) = 4; end % the minimum gap betwween overlapping bursts, in Hz. Default: a quarter of each frequency, with a minimum of 4
if ~isfield(pmtr,'max_span');    pmtr.max_span     = pmtr.foi(end)-pmtr.foi(1); end % the maximum allowed frequency span. Default: the whole spectrum
if ~isfield(pmtr,'wtpl_thr');    pmtr.wtpl_thr     = 0.5; end          % Threshold of WTPL to calculate duration
if ~isfield(pmtr,'min_dur');     pmtr.min_dur      = 1; end            % the minimum duration of a burst, in periods. Use 0.5 if the duration is calculated according to the high percentile
if ~isfield(pmtr,'ctrl_bursts'); pmtr.ctrl_bursts  = 'no'; end         % whether or not to plot individual bursts for control
if ~isfield(pmtr,'fs');         try pmtr.fs = raw.fsample;             % if fs (sapmling rate) not given, compute it  
                                catch pmtr.fs = round(1/nanmean(diff(raw.time{1}))); end; end 

disp('Calculating bursts');

% make the TFR
f_res = mode(diff(pmtr.foi)); % get the frequency resolution
foi = pmtr.foi;
nF = length(foi);
f1 = max([1, foi(1) - diff(foi([1,2]))]);
f2 = foi(end) + diff(foi([nF-1,nF]));
foi2 = [f1 pmtr.foi f2];
pmtr.min_f_gap = [pmtr.min_f_gap(1) pmtr.min_f_gap pmtr.min_f_gap(end)]; % add values to min_f_gap because we add a freq below and a freq above
disp(['Started calculating power at ' showCurrnetTime])
cfg             = [];
cfg.method      = 'wavelet'; %time-frequency analysis using the 'wavelet method' based on Morlet wavelets
cfg.output      = 'fourier'; %'pow'; % keep the full Fourier to support WTPL
cfg.foi         = foi2;
cfg.toi         = 'all';
cfg.width       = pmtr.wave_width;
cfg.pad         = 'nextpow2';
try TFR = ft_freqanalysis(cfg, raw);
catch
    warning('Having to remove next power of 2');
    cfg = rmfield(cfg,'pad');
    TFR = ft_freqanalysis(cfg, raw);
end
offline_pwr = squeeze(abs(TFR.fourierspctrm).^2);
% calculate WTPL
wfg                 = [];
wfg.verbose         = 1;
wfg.fs              = pmtr.fs;
wfg.nlag_post       = 1;
wfg.nlag_pre        = 1;
max_leg             = max([wfg.nlag_post, wfg.nlag_pre]);
min_start           = ceil(pmtr.fs/f1*max_leg)+1; % the minimum margin to keep from the ends
wfg.toi             = min_start:length(raw.time{1})-min_start;
WTPL = freqToWTPLwholeSession(wfg, TFR);
%
if ~isnan(pmtr.art_thresh)
    cfg                 = [];
    cfg.in_channel      = 1;
    cfg.fs              = pmtr.fs / pmtr.downsmpFact;
%     cfg.muscle          = pmtr.muscle;
    cfg.art_thresh      = pmtr.art_thresh;
    cfg.art_removal     = round(cfg.fs*0.5);
    pwr_clean = gross_artifact_removal(cfg, raw, offline_pwr);
    arts = isnan(gross_artifact_removal(cfg, raw, raw.trial{1}));
else 
    pwr_clean = offline_pwr;
    arts = [];
end

disp(['Started calculating burst peaks at ' showCurrnetTime])
t_BW = tic;
pwr_pctl_high = prctile (pwr_clean, pmtr.ptl_high, 2);
pwr_pctl_dur = prctile (pwr_clean, pmtr.ptl_dur, 2) + pmtr.SD_dur*std(pwr_clean,0,2,'omitnan');

pwr_masked = pwr_clean ./ pwr_pctl_high;
pwr_masked(pwr_masked < 1)=0;
pwr_masked(isnan(pwr_masked))=0;

pwr_for_dur = pwr_clean ./ pwr_pctl_dur;
pwr_for_dur(pwr_for_dur < 1)=0;
pwr_for_dur(isnan(pwr_for_dur))=0;

BW = imregionalmax(pwr_masked);
[i,j] = find(BW);
added = i==1 | i==length(foi2);
i(added) = [];
j(added) = [];
i(j>pmtr.lastSample)=[]; % remove bursts that happened after the task has ended
j(j>pmtr.lastSample)=[];

% translate to the original frequency frame
i = i-1;
BW ([1,end],:) = []; %one f above and one below were added just for peak reasons
pwr_clean([1,end],:) = [];
pwr_masked([1,end],:) = [];
pwr_for_dur([1,end],:) = [];
offline_pwr([1,end],:) = [];
pwr_pctl_dur([1,end],:) = [];
pwr_pctl_high([1,end],:) = [];
WTPL([1,end],:) = [];
thr = prctile(WTPL(:),75);
wtpl = WTPL;
wtpl(wtpl<thr)=0;

ii = i;% translate i into the coordinates of filtd
disp (['Calculating burst peaks took ' num2str(toc(t_BW)) ' seconds']);

% create array of filtered data, to detect the troughs and peaks
disp(['Started filtering at ' showCurrnetTime])
t_filter = tic;
raw_non_ft = raw.trial{1};
filtd = zeros(length(pmtr.foi),length(raw_non_ft));
for f = 1:length(pmtr.foi)
    try filtd(f,:) = ft_preproc_bandpassfilter(raw_non_ft, raw.fsample, pmtr.foi(f)+[-2 2]);
    catch bandpass(raw_non_ft,pmtr.foi(f)+[-1 1],raw.fsample); end
end
disp (['Filtering took ' num2str(toc(t_filter)) ' seconds']);

%%
disp(['Started burstDists structure at ' showCurrnetTime])
burstDists = zeros(length(j),21);
PeakSpan = cell(length(j),1);
PeakSamp = cell(length(j),1);
outTick = cell(length(j),2);% 1- timestapms of all troughs. 2= Peaks
% 1= maxima freq(Hz), 2= maxima timing(sample), 3= begin (sample), 4= end (sample), 5 = nearest trough
% 6 = duration (periods), 7 = duration (ms), 8 = power (rel pctl), 9 = power (absolute),
% 10 = number of troughs, 11 = number of peaks, 12= begin high percentile (sample), 13= end high percentile (sample)
% 14 = freq index (ignore), 15 = lowest SPAN frequency (Hz), 16 = highest SPAN frequency (Hz),
% 17 = lower PEAK freq (Hz), 18 = highest PEAK freq (Hz),
% 19 = begin (sample, wtpl), 20 = end(sample, wtpl), 21 = duration WTPL (periods)

too_short = false(length(j),1); % flag of keeping the burst
n = size(pwr_masked,2);
disp(['Total bursts: ' num2str(length(j))]);
t_trials = tic;
prev = fprintf(' ');
%
for x = 1:length(j)
    if ~mod(x,100); prev = dispRMVprev (['Finished ' num2str(x) ' bursts in ' num2str(toc(t_trials)) ' sec'],prev); end
    burstDists(x,1) = foi(i(x));% pmtr.foi(i(x)-1); %i(x)+pmtr.foi(1)-2;    % freq
    burstDists(x,14) = i(x);%i(x)+pmtr.foi(1)-2;    % freq index. to be removed
    burstDists(x,2) = j(x)*pmtr.downsmpFact; % timeStamp
    
    % get the beginning and end of the burst, according to the percentile used to detect the peak
    [t_before_high_pctl, t_after_high_pctl, ~, ~] = get_burst_durations(i,j,x,pmtr,ii,filtd,n,pwr_masked);
    % get the beginning and end of the burst, according to the percentile used to calculate duration
    [t_before, t_after, pks, trgh] = get_burst_durations(i,j,x,pmtr,ii,filtd,n,pwr_for_dur);
    
    burstDists(x,3) = t_before * pmtr.downsmpFact;
    burstDists(x,4) = t_after * pmtr.downsmpFact;
    burstDists(x,6) = (t_after - t_before) ./ (ceil(pmtr.fs./burstDists(x,1)));
    too_short(x) = burstDists(x,6)<pmtr.min_dur;
    
    if ~too_short(x)
        burstDists(x,7) = (t_after - t_before) / pmtr.fs*1000;
        burstDists(x,8) = pwr_masked(i(x), j(x));
        burstDists(x,9) = pwr_clean(i(x), j(x));
        tmp = filtd(ii(x), t_before-1:t_after+1);
        outTick{x,1} = trgh;
        outTick{x,2} = pks;
        burstDists(x,10) = length(trgh);
        burstDists(x,11) = length(pks);
        % finding the nearest trough
        wished_inds = ceil(pmtr.fs/burstDists(x,1)*1.5); % samples to look around the maxima, to make sure there is at least one trough
        pre_ind = max([1,j(x)-wished_inds]);
        post_ind = min([length(raw_non_ft),j(x)+wished_inds]);
        inds = pre_ind:post_ind;
        locs = [];
        val = [];
        for si = 2:(length(inds)-1) % this for loop is much faster than findpeaks
            if(filtd(ii(x), inds(si)) < filtd(ii(x), inds(si)-1) &&...
                    filtd(ii(x), inds(si)) < filtd(ii(x), inds(si)+1))
                locs = [locs, inds(si)]; val = [val, abs(filtd(ii(x), inds(si)))]; end
        end
        locs2 = abs(locs-burstDists(x,2));
        nearest_trough = locs(locs2==min(locs2));
        if length(nearest_trough)>1
            cand = val(locs2==min(locs2));
            [~,tmpi]=max(cand);
            nearest_trough = nearest_trough(tmpi);
        end
        if isempty(nearest_trough) % in the (very) rare condition that there is no trough, take the point of power maxima
            %                 figure; plot(filtd(ii(x), inds));
            nearest_trough = burstDists(x,2);
            %                 disp('Found a burst without troughs');
        end
        burstDists(x,5) = nearest_trough;
    end
    col = pwr_masked(:,burstDists(x,2));
    [freq_span,ind_above,ind_below] = getFreqSpan(col, i(x));
    burstDists(x,15) = foi(burstDists(x,14) - ind_below);
    burstDists(x,16) = foi(burstDists(x,14) + ind_above);
    % Peak frequency span
    B = burstDists(x,:);
    f_ind = B(14)-ind_below : B(14)+ind_above;
    if numel(f_ind)>1
        curr = pwr_masked(f_ind, B(3):B(4));
        [M,I] = max(curr); I(M==0) = [];
        burstDists(x,17) = foi(f_ind(min(I)));
        burstDists(x,18) = foi(f_ind(max(I)));
        peak_f_ind = f_ind(I);
        PeakSpan{x} = foi(peak_f_ind);
        PeakSamp{x} = B(3):B(4); PeakSamp{x}(M==0)=[];
    else
        burstDists(x,[17,18]) = foi(B(14)); % if there was only one frequency in the burst, the whole span is this frequency
        peak_f_ind = B(14);
        PeakSpan{x} = foi(peak_f_ind);
        PeakSamp{x} = B(1);
    end
    % Define the beginning and end of the burst according to all
    % frequencies that had a peak
    %         in_pwr = zeros(size(pwr_for_dur));
    %         in_pwr (B(14),:) = sum(pwr_for_dur(unique(peak_f_ind),:),1);
    %         get_burst_durations(i,j,x,pmtr,ii,filtd,n,in_pwr);
    %         in_pwr = sum(pwr_for_dur(unique(peak_f_ind),:),1);
    [t_before, t_after] = WTPL_burst_durations_v2(peak_f_ind, B, pwr_for_dur);
    burstDists(x,12) = t_before * pmtr.downsmpFact;
    burstDists(x,13) = t_after * pmtr.downsmpFact;
    burstDists(x,6) = (t_after - t_before) ./ (ceil(pmtr.fs./burstDists(x,1)));
    % WTPL-based duration
    % first aarguement is what freq indices to look at. Can be B(14) or unique(peak_f_ind)
    % B is the bursts distribution of this burst. wtpl should be normalized to threshold.
    %         in_wt = sum(wtpl(unique(peak_f_ind),:),1);
    [t_before, t_after] = WTPL_burst_durations_v2(peak_f_ind,B, wtpl);
    burstDists(x,19) = t_before;
    burstDists(x,20) = t_after;
    burstDists(x,21) = (t_after - t_before) ./ (ceil(pmtr.fs./burstDists(x,1)));
end
%%
bst = burstDists;
span = bst(:,16)- bst(:,15);
too_wide = span > pmtr.max_span;
span(too_wide)=[];
bst(too_short|too_wide,:)=[]; % remove bursts that were too short (time) or too wide (freq)
outTick(too_short|too_wide,:) = [];
PeakSpan(too_short|too_wide)=[];
%% Burst pruning
pwr = pwr_clean./pwr_pctl_high;
disp('.');
prev = fprintf(' ');
t_.bst=tic;
for bi = 1:length(bst)
    if ~mod(bi,2500); prev = dispRMVprev (['Finished pruning ' num2str(bi) ' bursts in ' num2str(toc(t_.bst)) ' sec'],prev); end
    %     if ~mod(bi,1e4); disp(['finished pruning ' num2str(bi) ' bursts in ' num2str(toc(t_.bst)) ' seconds']); end
    if bst(bi,1) > 0 % skip bursts that were already rejected
        cond1 = bst(:,3)>=bst(bi,3) & bst(:,3)<=bst(bi,4); % all bursts that begin after the beginning but before the end of the current
        cond2 = bst(:,4)<=bst(bi,4) & bst(:,4)>=bst(bi,3); % all bursts that finsh before the end but after the beginning of the current
        cond3 = bst(:,3)<=bst(bi,3) & bst(:,4)>=bst(bi,4); % all bursts that begin after and finish before the current
        cooc = cond1 | cond2 | cond3;
        cooc(bi)=0;
        cooc = find(cooc);
        for ci = 1:length(cooc)
            f1 = bst(bi,1);
            ind1 = bst(bi,14);
            f2 = bst(cooc(ci),1);
            ind2 = bst(cooc(ci),14);
            if abs(f1-f2) <= pmtr.min_f_gap(bst(bi,14)) % if the difference b/w the burst peakss is smaller than the defined minimum
                %                 pwr1 = pwr(ind1-pmtr.foi(1)+2, bst(bi,2)); % remove the burst with the lower energy.
                %                 pwr2 = pwr(ind2-pmtr.foi(1)+2, bst(cooc(ci),2)); % if they are equal, keep the first (which has the earlier peak).
                %                 % note that the power in the matrix is with 1 Hz below and
                %                 % 1 Hz above pmtr.foi.
                pwr1 = pwr(ind1, bst(bi,2)); % remove the burst with the lower energy.
                pwr2 = pwr(ind2, bst(cooc(ci),2)); % if they are equal, keep the first (which has the earlier peak).
                e1 = pwr1*(bst(bi,4)-bst(bi,3)); % estimation of energy as power x duration
                e2 = pwr2*(bst(cooc(ci),4)+bst(cooc(ci),3));
                if e2>e1
                    bst(cooc(ci),3) = min([bst(cooc(ci),3), bst(bi,3)]);
                    bst(cooc(ci),4) = max([bst(cooc(ci),4), bst(bi,4)]);
                    bst(bi,:) = 0; % reject the "bi" burst
                    break; % if the "original" was rejected no need to compare to the rest
                else
                    bst(bi,3) = min([bst(cooc(ci),3), bst(bi,3)]);
                    bst(bi,4) = max([bst(cooc(ci),4), bst(bi,4)]);
                    bst(cooc(ci),:) = 0; end
            end
        end
    end
end
outTick(~bst(:,1),:) = [];
PeakSpan(~bst(:,1),:) = [];
bst(~bst(:,1),:)=[];
disp('.');
disp (['Total bursts rejected: ' num2str(length(burstDists)- length(bst))]);
%%
Output.pmtr         = pmtr;
Output.burstDists   = bst;
Output.TroughPeak   = outTick;
Output.PeakSpan     = PeakSpan;
Output.arts         = arts;

Output.col_names    = {'maxima freq (Hz)','maxima timing (sample)',...
    'begin (sample)','end (sample)','nearest trough','duration (period)', ...
    'duration (ms)', 'power (rel pctl)', 'power (absolute)', 'number of troughs', 'number of peaks'...
    'begin high percentile (sample)', 'end high percentile (sample)', ...
    'ignore', 'low Span freq (Hz)', 'high Span freq (Hz)', 'low Peak freq (Hz)', 'high Peak freq (Hz)',...
    'begin (sample, wtpl)', 'end (sample, wtpl)', 'duration WTPL (periods)'}';


if pmtr.output_power 
    Output.pwr = pwr_clean./pwr_pctl_high;
    Output.filtd = filtd;
end
disp('.');
disp (['The call to Bursts_detection_WTPL took ' num2str(toc(t_func)) ,' seconds']);