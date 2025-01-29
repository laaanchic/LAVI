function PINK = computePinkLAVI(cfg,data)
% Generate pink noise to estimate the significance level of detected bands.
% Input:
% ======
%
% Data: N_channel x N_time array.
%
% cfg.Pink_reps : number of simulations created per channel. Default = 20.
% cfg.durs: the duration of each simulation. Default = duration of the data.

% Set the defaults
if ~isfield(cfg,'Pink_reps'); cfg.Pink_reps = 20; end % number of repetitions of pink noise simulations
if ~isfield(cfg,'foi'); cfg.foi = 10.^(0.5:0.025:1.65); end
if ~isfield(cfg,'fs'); cfg.fs = 1000; end
if ~isfield(cfg,'lag'); cfg.lag = 1.5; end
if ~isfield(cfg,'width'); cfg.width = 5; end
T = size(data,2)/cfg.fs;
if ~isfield(cfg,'durs'); cfg.durs = T; end % duration (in sec) of each simulation

if cfg.Pink_reps==0 || cfg.durs==0
    PINK = [];
else    
    Pink_reps   = cfg.Pink_reps;
    durs        = cfg.durs;
    foi         = cfg.foi;
    fs          = cfg.fs;
    lag         = cfg.lag;
    width       = cfg.width;
    pmtr        = cfg; % to keep a copy of the original cfg
    w =     hanning(fs*2); % 2-sec window
    
    N.time      = size(data,2);
    N.chan      = size(data,1);
    N.freq      = length(foi);
    
    PINK = nan(Pink_reps, N.freq, N.chan); % dimord = rep_freq_chan
    tmax = min([durs, N.time]);
    data = data(:,1:floor(tmax*fs)); % take shorter duration then the original if requested
    N.time      = size(data,2); % document the new number of samples
    
    t_.reps = tic;
    prev = fprintf('  ');
    for ch = 1:size(data,1)
        EEG = data(ch,:);
        EEG = EEG-mean(EEG); % demean
        %     sorted_values = sort(EEG); % use this when using the values of the original signal. Caution: sensitive to artifacts!
%         fourier_coeff = abs(ifft(EEG));
        [coefsIntoSurr, fito] = get_pink_iafft_coefs_pow (EEG,w,foi,fs);
        surrLAVI = nan(pmtr.Pink_reps,length(foi));
        rng = prctile(EEG,99.9)-prctile(EEG,0.01);%range(EEG); % use this  if using the option of random
        offset = prctile(EEG,0.01);%min(EEG); % values that are offsetted to the original EEG
        
        for ri = 1:pmtr.Pink_reps
            str = ['Running PINK ANALYSIS' ' Channel ' num2str(ch) '/' num2str(N.chan)...
                ' repeat ' num2str(ri) '/' num2str(pmtr.Pink_reps),...
                '. So far analysis took ' num2str(toc(t_.reps),4) ' seconds'];
            fprintf(repmat('\b',1,prev))
            prev = fprintf(str);
            
            %         pinkNoise = iaaft_loop_1d(coefsIntoSurr, sorted_values); % using the values of the original signal
            pinkNoise = iaaft_loop_1d(coefsIntoSurr, sort(rand(size(EEG)))); % using random signal and adjusting to range and offset
            %             pinkNoise = pinkNoise*rng + offset; % get the pink noise to the values of te EEG. Not super importnat for the LAVI analysis
            %         [pinkx,PINKfrq] = pwelch(pinkNoise,w,[],[],fs); % to document the powr of each pink simulation.2-sec window, whole session
            %         PINKpow(ri,:) = pinkx;%*pxx(i1)/pinkx(i1); % Uncomment if you want to document the power spectra of all pink simulations
            %         PINKnoise(ri,:) = pinkNoise; % Uncomment if you want to document all the pink noises
            
            for fi = 1:N.freq
                f = foi(fi);                
                spectrum = waveletLight(pinkNoise, fs, f, width);
                surrLAVI(ri,fi) = compute_lavi(spectrum,fs,f,lag);
            end
        end
        PINK(:,:,ch) = surrLAVI;
    end    
    disp('.');
end
end