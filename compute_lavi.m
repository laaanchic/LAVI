function [LAVI] = compute_lavi(spectrum,fs,f,lags)
% Generates the Lagged Angle Vector Index (LAVI) profile of one frequency
% over all channels in a session.
% syntax: [LAVI] = compute_lavi(spect,fs,f,lags)
%
% Mandatory input:
% ================
% spect: complex, N_chan x N_time, the Fourier spectrum of the signal
%        (as in the output of the function waveletLight.
% fs:   double, sampling frequency *of the transformed data* in Hz (that is,
%       the 1/time is seconds between each time bin).
% f:    frequency of interest, in Hz (should be the same as the freuqnecies of spect).
%       Better done one at a time, to avoid exceeding memory capacity.
%
% Optional input:
% ===============
% lags: the shift in time, measured in cycles, between the original signal
%       and the "lagged" copy. Default = 1.5.

if nargin<4; lags = 1.5; end
% nF = size(spectrum,2); % number of frequencies
nCH = size(spectrum,1);
nT = size(spectrum,3);
% if nF ~= length(f)
%     error('Number of frequencies in the data (dimension 2) and given frequencies do not match');
% end
LAVI = nan(nCH,1);
width = round(lags./f*fs); % lag size in samples
SIG0 = spectrum(:,1:end-width);
SIG1 = spectrum(:,width+1:end);
% remove NaNs
nanind = squeeze(isnan(SIG0(1,:)) | isnan(SIG1(1,:)));
SIG0(:,nanind) = [];
SIG1(:,nanind) = [];
A0 = abs(SIG0);
A1 = abs(SIG1);
for ch = 1:nCH
    sig0 = SIG0(ch,:);
    sig1 = SIG1(ch,:);
    a0 = A0(ch,:);
    a1 = A1(ch,:);
    LAVI(ch) = abs(sig0*sig1' /sqrt((a0*a0')*(a1*a1')));
end
