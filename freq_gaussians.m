clear

% specify wavelet parameters
peakf = 10;
fwhm  = 5.2;

% specify simulation details
npnts = 3000;
srate = 1000;

% vector of frequencies
hz = linspace(0,srate/2,npnts);

%% Morlet

% frequency-domain Gaussian (Cohen, 2019, Better Wavelets, Pg84 C.7)

s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-peakf;             % shifted frequencies
morlet_gaus = exp(-.5*(x/s).^2);    % gaussian

%% Morse

b = 27;
g = 3;

Wbg = (b/g)^(1/g); % maximum value at the peak frequency

% Normalisation factor (Olhede & Walden, 2002, P2666, k=0)
r=(2*b+1)/g;
A=sqrt(pi*g*2^r*exp(-gammaln(r)));

% frequency-domain Gaussian
y  = Wbg*hz;             % shifted frequencies
wa = 2*pi*y;
morse_gaus = sqrt(2)*A*(wa.^b).*(exp(-wa.^g));    % gaussian


%% plotting

time = (-floor(npnts/2):floor(npnts/2))/srate;

figure(1), clf
subplot(211)
plot(hz,morlet_gaus,'k','linew',2)
set(gca,'xlim',[0 peakf*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (gain)')

subplot(212), hold on
plot(hz,morse_gaus,'k','linew',2)
set(gca,'xlim',[0 peakf*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (gain)')