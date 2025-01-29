function spectrum = waveletLight(data, fsample, foi, width)
% A quick implementation of wavelet, based on ft_specest_wavelet
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

N.freq = length(foi);
N.time = size(data,2);
N.chan = size(data,1);
wltspctrm = cell(N.freq,1);
gwidth = 3; % that's the default used in fieldtrip
if isreal(data), signalFreq = fft(data);
else, signalFreq = data; end

for fi = 1:N.freq
    dt = 1/fsample;
    sf = foi(fi) / width;
    st = 1/(2*pi*sf);
    toi = -gwidth*st:dt:gwidth*st;
    A = 1/sqrt(st*sqrt(pi));
    tap = (A*exp(-toi.^2/(2*st^2)))';
    acttapnumsmp = size(tap,1);
    taplen(fi) = acttapnumsmp;
    ins = ceil(N.time./2) - floor(acttapnumsmp./2);
    prezer = zeros(ins,1);
    pstzer = zeros(N.time - ((ins-1) + acttapnumsmp)-1,1);
    
    % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
    ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./fsample) .* foi(fi));
    
    % create wavelet and fft it
    wavelet = complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer));
    wltspctrm{fi} = complex(zeros(1,N.time));
    wltspctrm{fi} = fft(wavelet,[],1)';
end

timeboi = 1:N.time;
spectrum = complex(nan(N.chan,N.freq,N.time),nan(N.chan,N.freq,N.time));
for fi = 1:N.freq
    % compute indices that will be used to extracted the requested fft output
    nsamplefreqoi    = taplen(fi);
    reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi < (N.time - (nsamplefreqoi ./2))));
    reqtimeboi       = timeboi(reqtimeboiind);
    if ~isempty(reqtimeboi)   
        dum = fftshift(ifft(signalFreq .* repmat(wltspctrm{fi},[N.chan 1]), [], 2),2);
%         dum = fftshift(ifft(signalFreq .* wltspctrm{fi}, [], 2),2);        
        dum = dum .* sqrt(2 ./ fsample);
        spectrum(:,fi,reqtimeboiind) = dum(:,reqtimeboi);
    end
end