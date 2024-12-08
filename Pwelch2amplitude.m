function amp = Pwelch2amplitude(pxx,f,w)
% Syntax: amp = Pwelch_to_amplitude(pxx,f,w)
% Input:
% pxx and f: the power and frequency output of pwelch
% w:         the window (input) to pwelch
% Translates the PSD obtained from pwelch into amplitude.
% Based on the document
% "How to use the FFT and Matlab’s pwelch function for signal and noise simulations and measurements"
% google: amplitude estimation using PWelch
% http://schmid-werren.ch/hanspeter/publications/ 
fbin = f(2)-f(1);
CG = sum(w)/length(w); 
NG = sum(w.^2)/length(w);
amp = sqrt(pxx * (NG*fbin)/(CG^2) * 2);