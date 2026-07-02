function [LAVI,cfg] = Prepare_LAVI(cfg,data)
% Takes raw data and transfers it to the frequency domain, as preparation
% for the compute_lavi function
% syntax: [LAVI,cfg] = Prepare_LAVI(cfg,data)
% configuration fields:
% cfg.foi =     frequencies of interest. default: 10.^(0.5:0.025:1.65)
% cfg.fs =      sampling frequency. Default: 1000 Hz.
% cfg.lag =     the time delay between the data and the copy of itself, in 
%               cycles. Default: 1.5.
% cfg.width =   the width, in cycles, of the wavelet. Default: 5.
% cfg.verbose = whether display messages on screen. Default: 1 (yes).

tic;
if ~isfield(cfg,'foi'); cfg.foi = 10.^(0.5:0.025:1.65); end
if ~isfield(cfg,'fs'); cfg.fs = 1000; end
if ~isfield(cfg,'lag'); cfg.lag = 1.5; end
if ~isfield(cfg,'width'); cfg.width = 5; end
if ~isfield(cfg,'verbose'); cfg.verbose = 1; end

st = dbstack;
funcname = st.name; clear st;

foi         = cfg.foi;
fs          = cfg.fs;
lag         = cfg.lag;
width       = cfg.width;

N.time      = size(data,2);
N.chan      = size(data,1);
N.freq      = length(foi);
%% Calculate wavelet and LAVI per frequency
LAVI = zeros(N.chan, N.freq);
if cfg.verbose, prev = fprintf(' '); end
for fi = 1:N.freq
    f = foi(fi);
    if cfg.verbose
        str = ['Running LAVI frequency ' num2str(fi) '/' num2str(N.freq) ' (' num2str(f,3) ' Hz)'];
        fprintf(repmat('\b',1,prev))
        prev = fprintf(str);
    end
    for ch = 1:N.chan
        if any(isnan(data(ch,:)))
            spectrum = squeeze(tfrLight(data(ch,:), fs, f, width, 0));            
        else
            spectrum = squeeze(waveletLight(data(ch,:), fs, f, width));
        end
        if iscolumn(spectrum), spectrum=spectrum'; end
        lavi = compute_lavi(spectrum,fs,f,lag);
        LAVI(ch,fi) = lavi;
    end
end

if cfg.verbose
    disp('.');
    disp (['The call to ' funcname ' took ' num2str(toc,3) ' seconds']);
end
end