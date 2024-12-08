function [BORDERS,varNames,SIGVECT] = ABBA(LAVI,foi,alpha_range,SIGLIM,perFreq)
% Finds the bands, borders, and significance based on the LAVI profile.
%
% Syntax: [BORDERS,varNames,SIGVECT] = ABBA(LAVI,foi,alpha_range,sigLim,perFreq)
%
% Input:
% ======
% LAVI: the output of the function compute_lavi
%
% foi: the frequencies of LAVI.
%
% alpha_range: the frequency range in which we expact to find alpha. The
%      band of the peak in this band will be assigned the index 0, bands
%      with lower frequency will be assigned with negative indices, and
%      above with positive indices.
%
% SIGLIM: N_chan x N_freq x 2 array, with the lower and upper values of 
%      significance level, which can be obtained, e.g., by a distribution 
%      of pink noise. Default: the median of the LAVI profile.
%
% perFreq: bool, how to define the significance level: 
%      1= per frequency 
%      0= min/ max over all frequncies 
%
% Output:
% =======
% BORDERS:  1 x N_channels cell array. Each array contains an N_bands x 11
%           matrix, with the following columns (also documented in varNames):
%           1. Index of band beginning.
%           2. Index of band end.
%           3. Index of the band peak/ trough.
%           4. Frequency at the beginning of the band.
%           5. Frequency at the end of the band.
%           6. Frequency at the peak of the band.
%           7. The LAVI value at the the peak.
%           8. The LAVI value relative to the median at the peak.
%           9. The direaction of the band (1= peak, -1= trough).
%           10. Index of the band relative to alpha. Negative: frequencies
%               lower than alpha. Positive: frequencies higher than alpha.
%               Useful to assign names to bands.
%           11. Boolean, whether the peak/trough of the band is significant.
%
% SIGVECT: 1 x N_channels cell array. Each array contains a 1 X N_freqs
%          array. Frequencies with significantly high LAVI are assigned
%          with positive numbers. Frequencies with significantly low LAVI 
%          are assigned with negative numbers. Non-significant frequencies
%          are assigned with 0.

if nargin<3 || isempty(alpha_range), alpha_range = [6 14]; end

varNames = {'BegI','EndI','PeakI','BegF','EndF','PeakF','PeakLAVI','PeakRel','Dir','Rel_alpha','Sig'}';
N.chan = size(LAVI,1);
N.freq = size(LAVI,2);
BORDERS = cell(1, N.chan); % the borders of each channel will be saved in a cell (because there can be different numbers of bands in different channels)
SIGVECT = cell(1, N.chan);
if size(alpha_range,1)==1, alpha_range = repmat(alpha_range,N.chan,1); end
for ch = 1:N.chan
    lavi = LAVI(ch,:);    
    if nargin<4 || isempty(SIGLIM)
        sigLim = repmat([0;0]+nanmedian(lavi),1,N.freq);
    elseif ndims(SIGLIM)==3, sigLim = squeeze(SIGLIM(ch,:,:));
%     elseif size(SIGLIM,1)>1, sigLim = SIGLIM(ch,:);
    elseif numel(SIGLIM)>2, sigLim = SIGLIM;
    else % there is an option to provide only minimum and maximum values
        if isrow(SIGLIM), SIGLIM=SIGLIM'; end
        sigLim = repmat(SIGLIM,1,N.freq);
    end
    if size(sigLim,1)~=2,sigLim=sigLim';end % ensure the min/max is the first dimension
    if size(sigLim,1)~=2, error('Wrong definition of SIGLIM'); end
    if nargin<4 || ~perFreq % choose whether to use siginificance level per frequency, or as the minimum/ maximum over all frequencies
        sigLim(1,:) = min(sigLim(1,:));
        sigLim(2,:) = max(sigLim(2,:));
    end
    % significance vector
    sigVect = zeros(size(foi));
    sigVect(lavi>sigLim(2,:)) = 0.5;
    sigVect(lavi<sigLim(1,:)) = -0.5;
    
    % find band limits as crossings of the reference. 
    % Deal with (rare) cases of points equal exactly to the reference
    ref = nanmedian(lavi);
    reref = lavi'-ref;
    reref = dealWithZeros(reref);
    siman = sign(reref);
    flipp = [diff(siman);0]; % the last freq in each band is non zero   
    flipp(isnan(flipp)) = 0;
    borders = find(flipp); % the borders are where the sign has changed    
    borders = [[1;borders+1], [borders;length(foi)]];
    if borders(end,1)>length(foi); borders(end,:)=[]; end    
    borders(:,9) = siman(borders(:,2));    
        
    for bi = 1:size(borders,1)
        inds = borders(bi,1):borders(bi,2);
        [~,ind] = max(abs(reref(inds)));
        ind = ind+inds(1)-1;
        borders(bi,3) = ind;
        borders(bi,7) = lavi(ind);
        borders(bi,8) = reref(ind);
    end
    borders(:,4:6) = round(foi(borders(:,1:3)),1);
    % defining alpha as band 0    
%     w_inds = find(borders(:,4)<=alpha_range(ch,2) & borders(:,5)>=alpha_range(ch,1)); % bands with width falling within defined alpha range
    w_inds = find(borders(:,6)<=alpha_range(ch,2) & borders(:,6)>=alpha_range(ch,1) & borders(:,9)>0); % bands with peaks falling within defined alpha range
    if ~isempty(w_inds)
        [~,I] = max(borders(w_inds,7));
        a_ind = w_inds(I);
        borders(:,10) = (1:size(borders,1))'-a_ind;
    else
%         [~,forcedAlphaInd] = max(lavi.*(foi<=alpha_range(ch,2)&foi>=alpha_range(ch,1)));
        borders(:,[4:6,10]) = nan;
    end
    borders(:,11) = (sigVect(borders(:,3))'.*borders(:,9))==0.5; % detemines if the peak of the band is signiifcant
    nsigb = find(~borders(:,11));
    for ni = 1:length(nsigb)
        bi = nsigb(ni);
        sigVect(borders(bi,1):borders(bi,2)) = 0;
    end
    BORDERS{ch} = borders;
    SIGVECT{ch} = sigVect;
end

function output = dealWithZeros(input)
    output = input;
    zs = find(input==0);
    if ~isempty(zs)
        for i = 1:length(zs)
            if zs(i)==1, output(1)=output(2)/10; 
            else, output(zs(i)) = output(zs(i)-1)/10; % keep the sign of the previous sample, but get closer to 0
            end
        end
    end
end
end