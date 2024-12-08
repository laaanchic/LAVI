function [coefsIntoSurr, amp] = get_pink_iafft_coefs_random_ap (N,fs,foi,a,b)
% Syntax: [coefsIntoSurr, amp] = get_pink_iafft_coefs_random_ap (N,fs,ap)
% Generates the coefficients to be used to simulate pink noise using iafft
w = hanning(fs*2);
EEG = rand(1,N);
[~,pff] = pwelch(EEG,w,[],[],fs);
ff = getFrequenciesOfFFT(fs, N);
posf = ff(ff>0);
pow = a.*pff.^b;
amp = Pwelch2amplitude(pow,pff,w);
amp = amp*N/2; % amplitude to coefficients

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