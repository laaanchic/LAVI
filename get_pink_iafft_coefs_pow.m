function [coefsIntoSurr, fito, amp] = get_pink_iafft_coefs_pow (EEG,w,foi,fs)
% Syntax: [coefsIntoSurr, fito, amp] = get_pink_iafft_coefs_pow (EEG,w,foi,fs)
% Generates the coefficients to be used to simulate pink noise using iafft
% fito is the fit object
N = length(EEG);
[pxx,pff] = pwelch(EEG',w,[],[],fs);
ff = getFrequenciesOfFFT(fs, length(EEG));
posf = ff(ff>0);
amp = Pwelch2amplitude(pxx,pff,w);
amp = amp*N/2; % amplitude to coefficients
% amp = pxx;

fitFs = pff>=foi(1) & pff<=foi(end);
fitx = double(pff(fitFs));
fity = double(amp(fitFs));
fito = fit(fitx, fity, 'power1');
a=fito.a; b=fito.b; %c=fito.c;
ap = (a.*posf.^b); 
switch mod(N,2)
    case 1
        coefsIntoSurr = [0, ap, flip(ap)];
    case 0
        coefsIntoSurr = [0, ap, flip(ap(1:end-1))];
end