function [spectrum, foi, toi] = tfrLight(data, fsample, foi, width, verbose)
% (dat, time, varargin) input in the original function
% (data, fsample, foi, width) input in my waveletLight
% A quick implementation of tfr (convoultion of wavelet in the time domain),
% based on ft_specest_tfr
%
% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% set the default for screen output
try ft_warning off; end
if nargin<5, verbose=1; end

% Set n's
[nchan,ndatsample] = size(data);
% de-mean
data = data-nanmean(data);
% define time of interest. Always take the whole session
time = (1:ndatsample)/fsample;
toi = time;
gwidth = 3; % that's the default if ft

% set total time-length of data
dattime = ndatsample / fsample; % total time in seconds of input data
pad = 2^nextpow2(dattime);
% Zero padding
if round(pad * fsample) < ndatsample
    error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
    pad = dattime;
end
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data
prepad     = zeros(1,floor(((pad - dattime) * fsample) ./ 2));
postpad    = zeros(1,ceil(((pad - dattime) * fsample) ./ 2));

% Set freqboi and freqoi
freqoiinput = foi;
if isnumeric(foi) % if input is a vector
    freqboi   = round(foi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
    freqboi   = unique(freqboi);
    foi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(foi,'all') % if input was 'all'
    freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
    freqboi    = freqboilim(1):1:freqboilim(2);
    foi     = (freqboi-1) ./ endtime;
end
% check for freqoi = 0 and remove it, there is no wavelet for freqoi = 0
if foi(1)==0
    foi(1)  = [];
    freqboi(1) = [];
end
nfreqboi = length(freqboi);
nfreqoi  = length(foi);

% % throw a warning if input freqoi is different from output freqoi
% if isnumeric(freqoiinput)
%     % check whether padding is appropriate for the requested frequency resolution
%     rayl = 1/endtime;
%     if any(rem(freqoiinput,rayl)) % not always the case when they mismatch
%         warning('padding not sufficient for requested frequency resolution, for more information please see the FAQs on www.ru.nl/neuroimaging/fieldtrip');
%     end
%     if numel(freqoiinput) ~= numel(foi) % freqoi will not contain double frequency bins when requested
%         warning('output frequencies are different from input frequencies, multiples of the same bin were requested but not given');
%     else
%         if any(abs(freqoiinput-foi) >= eps*1e6)
%             warning('output frequencies are different from input frequencies');
%         end
%     end
% end


% Set timeboi and timeoi
timeoiinput = toi;
offset = round(time(1)*fsample);
if isnumeric(toi) % if input is a vector
    toi   = unique(round(toi .* fsample) ./ fsample);
    timeboi  = round(toi .* fsample - offset) + 1;
    ntimeboi = length(timeboi);
elseif strcmp(toi,'all') % if input was 'all'
    timeboi  = 1:length(time);
    ntimeboi = length(timeboi);
    toi   = time;
end

% % throw a warning if input timeoi is different from output timeoi
% if isnumeric(timeoiinput)
%     if numel(timeoiinput) ~= numel(toi) % timeoi will not contain double time-bins when requested
%         warning('output time-bins are different from input time-bins, multiples of the same bin were requested but not given');
%     else
%         if any(abs(timeoiinput-toi) >= eps*1e6)
%             warning('output time-bins are different from input time-bins');
%         end
%     end
% end


% Creating wavelets
% expand width to array if constant width
if numel(width) == 1
    width = ones(1,nfreqoi) * width;
end
wavelet = cell(nfreqoi,1);
for ifreqoi = 1:nfreqoi
    dt = 1/fsample;
    sf = foi(ifreqoi) / width(ifreqoi);
    st = 1/(2*pi*sf);
    toi2 = -gwidth*st:dt:gwidth*st;
    A = 1/sqrt(st*sqrt(pi));
    tap = (A*exp(-toi2.^2/(2*st^2)))';
    acttapnumsmp = size(tap,1);
    taplen(ifreqoi) = acttapnumsmp;
    ins = ceil(endnsample./2) - floor(acttapnumsmp./2);
    %prezer = zeros(ins,1);
    %pstzer = zeros(endnsample - ((ins-1) + acttapnumsmp)-1,1);
    
    % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
    ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./fsample) .* foi(ifreqoi));
    
    % create wavelet and fft it
    %wavelet{ifreqoi} = complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer));
    wavelet{ifreqoi} = complex(vertcat(tap.*cos(ind)), vertcat(tap.*sin(ind)));
end


% compute spectrum by convolving the wavelets with the data
spectrum = complex(nan(nchan,nfreqoi,ntimeboi),nan(nchan,nfreqoi,ntimeboi));
for ifreqoi = 1:nfreqoi
    if verbose
        str = sprintf('frequency %d (%.2f Hz)', ifreqoi,foi(ifreqoi));
        [st, cws] = dbstack;
        fprintf([str, '\n']);
    end
    
    
    
    % compute indices that will be used to extracted the requested output (this keeps nans when the wavelet is not fully immersed in the data)
    nsamplefreqoi    = taplen(ifreqoi);
    reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi < (ndatsample - (nsamplefreqoi ./2))));
    reqtimeboi       = timeboi(reqtimeboiind);
    
    % do convolution, if there are reqtimeboi's that have data
    if ~isempty(reqtimeboi)
        dum = complex(zeros(nchan,numel(reqtimeboi)));
        for ichan = 1:nchan
            dumconv = conv(data(ichan,:),  wavelet{ifreqoi}, 'same');
            dum(ichan,:) = dumconv(reqtimeboi); % keeping nans nans when the wavelet is not fully immersed in the data
        end
        spectrum(:,ifreqoi,reqtimeboiind) = dum;
    end
end
try ft_warning on; end
end

