function [PXX, f] = pwelchNaN(X, window, noverlap, nfft,fs)
% Welch's Power Spectral Density (PSD) estimation without using pwelch()
%
% Inputs:
%   x        - Input signal
%   fs       - Sampling frequency (Hz)
%   window   - Window function (vector)
%   noverlap - Number of overlapping samples
%   nfft     - Number of FFT points
%
% Outputs:
%   Pxx - Power spectral density estimate
%   f   - Frequency vector

L = length(window);  % Window length
if nargin<3 || isempty(noverlap), noverlap = L/2; end
if nargin<4 || isempty(nfft), nfft = 2^nextpow2(L); end
step = L - noverlap; % Step size for segmenting
n_chan = size(X,2);
half_idx = floor(nfft/2) + 1;
PXX = zeros(half_idx, n_chan);

for chi = 1:n_chan
    x = X(:,chi);
    num_segments = floor((length(x) - noverlap) / step); % Number of segments
    
    % Initialize PSD accumulator    
    Pxx = zeros(nfft, num_segments);
    % Power correction factor due to windowing
    U = sum(window.^2) / L;
    rmv_seg = [];
    % Loop over segments
    for i = 1:num_segments
        idx_start = (i-1) * step + 1;
        idx_end = idx_start + L - 1;
        
        % Extract segment and apply window function
        if idx_end <= length(x)
            segment = x(idx_start:idx_end) .* window;
        else
            break;
        end
        if any(isnan(segment))
            rmv_seg = [rmv_seg,i]; % mark segments to be removed
        else % Compute FFT of the segment
            Xf = fft(segment, nfft);
            
            % Compute power spectrum of the segment
            Pxx(:,i) = (abs(Xf).^2) / (L * fs);
            %             Pxx = Pxx + (abs(X).^2) / (L * fs);
        end
    end
    
    % Average over all segments and normalize by window power
    Pxx(:,rmv_seg) = [];
    num_segments = num_segments - length(rmv_seg);
    Pxx = sum(Pxx,2) / (num_segments * U);
    
    % Generate frequency axis
    f = (0:nfft-1) * (fs / nfft);
    
    % Keep only the one-sided spectrum (for real signals)    
    Pxx = Pxx(1:half_idx);
    % Scale non-DC and non-Nyquist frequencies by 2
    Pxx(2:end-1) = 2 * Pxx(2:end-1);
    f = f(1:half_idx);
    PXX(:,chi) = Pxx;
end
end
