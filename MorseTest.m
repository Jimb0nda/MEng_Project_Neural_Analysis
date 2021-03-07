%% Setup Parameters
srate = 1000;           % in Hz
wavelt = 0:1/srate:35;  % best practice to have 0 in the center of the wavelet

%Wavelet Parameters

% Frequency parameters
min_freq = 10; % Hz
max_freq = 35; % Hz
num_freq = 40; % count
frex = linspace(min_freq,max_freq,num_freq);

unit_step = wavelt>=0;
a = 0.501;
gamma = 3;
beta = 27;

Wbg = (beta/gamma)^(1/gamma); % maximum value at the peak frequency

%% Morse Wavelet

awt = unit_step .* a .* ((2*pi.*wavelt).^beta) .* (exp(-(2*pi.*wavelt).^gamma));

x = ifft(awt);
%% Plotting

figure(3), clf
plot(wavelt,imag(x));
hold on
title("Beta = " + beta + " & Gamma = " + gamma)
hold off
