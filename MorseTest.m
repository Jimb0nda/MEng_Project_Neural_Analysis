clear

%% Setup Parameters
srate = 1000;           % in Hz

%Wavelet Parameters
a = 0.501;
gamma = 3;
beta = 27;

Wbg = (beta/gamma)^(1/gamma); % maximum value at the peak frequency

% Frequency parameters
min_freq = 0; % Hz
max_freq = Wbg; % Hz
num_freq = 3000; % count
frex = linspace(min_freq,max_freq,num_freq);

%% Morse Wavelet

awt =  a .* ((2*pi.*frex).^beta) .* (exp(-(2*pi.*frex).^gamma));

x = ifft(awt);
%% Plotting

figure(4), clf
subplot(121);
plot(frex,real(x));
hold on
plot(frex,imag(x));
title("Beta = " + beta + " & Gamma = " + gamma)
hold off

subplot(122)
plot(frex,awt);
title("Beta = " + beta + " & Gamma = " + gamma)

