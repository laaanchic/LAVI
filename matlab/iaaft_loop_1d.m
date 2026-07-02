function [y, reps, errorAmplitude, errorSpec] = iaaft_loop_1d(fourierCoeff, sortedValues, maxReps)
% INPUT: 
% fourierCoeff:   The 1 dimensional Fourier coefficients that describe the
%                 structure and implicitely pass the size of the matrix
% sortedValues:   A vector with all the wanted amplitudes (e.g. LWC of LWP
%                 values) sorted in acending order.
% maxReps:        An integer, maximal number of repetitions before exiting
%                 the function if parametrical cinditions are not met
% OUTPUT:
% y:              The 1D IAAFT surrogate time series
% errorAmplitude: The amount of addaption that was made in the last
%                 amplitude addaption relative to the total standard deviation.
% errorSpec:      The amont of addaption that was made in the last fourier
%                 coefficient addaption relative to the total standard deviation
% When using this script, please credit the original contribution:
% V. Venema (2023). Surrogate time series and fields 
% (https://www.mathworks.com/matlabcentral/fileexchange/4783-surrogate-time-series-and-fields), 
% MATLAB Central File Exchange. Retrieved January 17, 2023. 
%
% License:
% Copyright (c) 2003, V. Venema
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Settings
errorThresshold = 2e-4; %
timeThresshold  = Inf;  % Time in seconds or Inf to remove this condition
speedThresshold = 1e-6; % Minimal convergence speed in the maximum error.
verbose = 0; % comments on screen
makePlots = 0; % Best used together with debugging.

% Initialse function
% noValues = size(fourierCoeff);
errorAmplitude = 1;
errorSpec = 1;
if nargin<3, maxReps = []; end
if isempty(maxReps), maxReps = 1000; end % maximal number of repetitions
oldTotalError = 100;
speed = 1;
standardDeviation = std(sortedValues);
t = cputime;

% The method starts with a randomized uncorrelated time series y with the pdf of
% sorted_values
[~,index]=sort(rand(size(sortedValues)));
y(index) = sortedValues;
reps = 0;

% Main intative loop
while ( (errorAmplitude > errorThresshold || errorSpec > errorThresshold) && (cputime-t < timeThresshold) && (speed > speedThresshold) && reps < maxReps)
    reps = reps + 1; % repetition counter
    % addapt the power spectrum
    oldSurrogate = y;    
    x=ifft(y);
    phase = angle(x);
    x = fourierCoeff .* exp(1i*phase);
    y = fft(x);
    difference=mean(mean(abs(real(y)-real(oldSurrogate))));
    errorSpec = difference/standardDeviation;
    if ( verbose ), errorSpec, end
    
    if (makePlots)
        plot(real(y))
        title('Surrogate after spectal adaptation')
        axis tight
        pause(0.01)
    end
        
    % addept the amplitude distribution
    oldSurrogate = y;
    [~,index]=sort(real(y));
    y(index)=sortedValues;
    difference=mean(mean(abs(real(y)-real(oldSurrogate))));
    errorAmplitude = difference/standardDeviation;
    if ( verbose ), errorAmplitude, end
    
    if (makePlots)
        plot(real(y))
        title('Surrogate after amplitude adaptation')
        axis tight
        pause(0.01)
    end
    totalError = errorSpec + errorAmplitude;
    speed = abs((oldTotalError - totalError) / totalError);
    if ( verbose ), totalError, speed, end
    oldTotalError = totalError;
end
% errorSpec
% errorAmplitude
y = real(y);
