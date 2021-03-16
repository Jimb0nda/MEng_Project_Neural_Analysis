clear

%% Setup Parameters
srate = 1000;           % in Hz
npnts = 3001;

% vector of frequencies
hz = linspace(0,35,npnts);

b = 9;
g = 3;

Wbg = (b/g)^(1/g); % maximum value at the peak frequency

% Normalisation factor (Olhede & Walden, 2002, P2666, k=0)
r=(2*b+1)/g;
A=sqrt(pi*g*2^r*exp(-gammaln(r)));

y  = Wbg*hz;             % shifted frequencies
wa = 2*pi*y;

% Normalisation factor for unit energy at each scale
scale_fac=sqrt(srate*Wbg);

%% Morse Wavelet

awt = scale_fac*sqrt(2)*A*(wa.^b).*(exp(-wa.^g));
awt = awt ./ max(awt);

x = fftshift(ifft(awt));
%% Plotting

time = (-floor(npnts/2):floor(npnts/2))/srate;

figure(4), clf
subplot(121);
plot(time,real(x));
hold on
plot(time,imag(x));
title("Time Domain Morse Wavelet: β = " + b + " & γ = " + g)
hold off

subplot(122)
plot(hz,awt);
title("Freq Domain Morse Wavelet: β = " + b + " & γ = " + g)

