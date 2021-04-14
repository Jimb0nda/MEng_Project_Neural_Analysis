% parameters
srate = 1000;         % in hz
time  = -1:1/srate:1; % best practice is to have time=0 at the center of the wavelet
frex  = 2*pi;         % frequency of wavelet, in Hz

% create sine wave (actually cosine, just to make it nice and symmetric)
sine_wave = cos( 2*pi*frex.*time );

% create Gaussian window
fwhm = .3; % width of the Gaussian in seconds
fwhm2 = .7; % width of the Gaussian in seconds
gaus_win = exp( (-4*log(2)*time.^2) / (fwhm^2) );

% create wavelet and get its FFT
cmw = exp(2*1i*pi*4.*time) .* exp( (-4*log(2)*time.^2) / (fwhm^2) ); % Morlet Wavelet
% create wavelet and get its FFT
cmw2 = exp(2*1i*pi*4.*time) .* exp( (-4*log(2)*time.^2) / (fwhm2^2) ); % Morlet Wavelet

pnts = length(time);

mwX = abs(fft( cmw )/pnts);
mwX2 = abs(fft( cmw2 )/pnts);
hz  = linspace(0,srate,pnts);


figure(1), clf

subplot(221), hold on
plot(time,real(cmw),'k')
xlabel('Time (s)'), ylabel('Amplitude')
title('Morlet wavelet with 3 cycles')

subplot(222), hold on
plot(time,real(cmw2),'r')
xlabel('Time (s)'), ylabel('Amplitude')
title('Morlet wavelet with 7 cycles')

subplot(2,2,[3,4]), hold on
plot(hz,mwX,'k')
plot(hz,mwX2,'r')
set(gca,'xlim',[0 frex*2])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
legend({'3 Cycles';'7 cycles'})
title('Complex Morlet wavelet in the frequency domain')