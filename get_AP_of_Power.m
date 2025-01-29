function [a,b] = get_AP_of_Power(data,fs,flim)
% Syntax: [a,b] = get_AP_of_Power(data,fs,flim)
% use this function to find the offset (a) and slope (b) of the aperiodic
% component of the POWER of data.
% These values can be used to generate pink noise with get_pink_iafft_coefs_random_ap
if nargin<3, flim = [5 40]; end
fs = round(fs);
winsize = min([length(data), 2*fs]);
w = hanning(winsize);
[pxx,pff] = pwelch(data,w,[],[],fs);
indf = pff>=flim(1) & pff<=flim(end);
% %%
% figure(7523);clf; set(gcf,'position',[490 100 560 420])
% plot(pff(indf),pxx(indf)); %xlim([2 40])
% %%
% fito = fit(pff(indf), pxx(indf), 'power1');
% a = fito.a; 
% b = fito.b; 
x = pff(indf);
y = pxx(indf);
f = @(a,b,x) a.*x.^b;
fito = @(params) norm(f(params(1), params(2),x)-y);
sol = fminsearch(fito, [1,1]);
a = sol(1);
b = sol(2);
