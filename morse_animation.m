clear

% Setup Parameters
srate = 1000;           % in Hz
npnts = 3000;

freq = linspace(0,srate/2,npnts);

% Frequency parameters for Morse
min_freq = 4; % Hz
max_freq = 35; % Hz
num_freq = 40; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% movie time parameter!
refresh_speed = .6; % seconds

b = 9;
g = 3;

Wbg = (b/g)^(1/g); % maximum value at the peak frequency

% Normalisation factor (Olhede & Walden, 2002, P2666, k=0)
r=(2*b+1)/g;
A=sqrt(pi*g*2^r*exp(-gammaln(r)));

for i = 1:num_freq

    a = Wbg/frex(i);
    y  = a*freq;        
    wa = 2*pi*y;

    % Normalisation factor for unit energy at each scale
    scale_fac=sqrt(srate*Wbg);

    % Generate normalised, scaled wavelet (Olhede & Walden, 2002, P2666, (10), k=0)
    awt = (scale_fac*sqrt(2)*A*(wa.^b).*(exp(-wa.^g)))'; 

    x = fftshift(ifft(awt));
    
    time = (-floor(npnts/2):floor(npnts/2))/srate;

    figure(4), clf
    subplot(121);
    plot(real(x));
    hold on
    plot(imag(x));
    title("Time Domain Morse Wavelet: β = " + b + " & γ = " + g)
    hold off

    subplot(122)
    plot(freq,awt);
    title("Freq Domain Morse Wavelet: β = " + b + " & γ = " + g)

    
    pause(refresh_speed)
end


