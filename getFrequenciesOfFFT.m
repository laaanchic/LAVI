function f = getFrequenciesOfFFT(fs, N)
% f = getFrequenciesOfFFT(fs, N)
% gets the frequencies after FFT, according to the Nyquist frequency and
% the number of samples. Takes care of both odd and even number of samples
f = (0:N-1)*fs/N;
switch mod(N,2)
    case 1
        f((N+1)/2+1:end) = f((N+1)/2+1:end)-fs;
    case 0
        f(N/2+2:end) = f(N/2+2:end)-fs;
end