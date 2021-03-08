clear

%% Setup Parameters
srate = 1000;           % in Hz

%Wavelet Parameters
a = 0.501;
gamma = 3;
beta = 0.333;

Wbg = (beta/gamma)^(1/gamma); % maximum value at the peak frequency

% Frequency parameters
min_freq = 0; % Hz
max_freq = Wbg; % Hz
num_freq = 100; % count
frex = linspace(min_freq,max_freq,num_freq);

%% Morse Wavelet

awt =  a .* ((2*pi.*frex).^beta) .* (exp(-(2*pi.*frex).^gamma));
awt = awt ./ max(awt)

x = fftshift(ifft(awt));
%% Plotting

figure(4), clf
subplot(121);
plot(frex,real(x));
hold on
plot(frex,imag(x));
title("Time Domain Morse Wavelet: β = " + beta + " & γ = " + gamma)
hold off

subplot(122)
plot(frex,awt);
title("Freq Domain Morse Wavelet: β = " + beta + " & γ = " + gamma)

