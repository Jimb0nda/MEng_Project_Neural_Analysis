clear

%% Setup Parameters
srate = 1000;           % in Hz
npnts = 20000;

hz = linspace(0,srate,npnts);

% Frequency parameters for Morse
min_freq = 4; % Hz
max_freq = 35; % Hz
num_freq = 40; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% Define complete frequency range
freq_all=[0:npnts/2,-npnts/2+1:-1]'*(srate/npnts);
% Apply heaviside function, only positive frequencies, first (T/2)+1 pts.
freq=freq_all(find(freq_all>=0));

g = 3;
b = 3;

Wbg = (b/g)^(1/g); % maximum value at the peak frequency

% Normalisation factor (Olhede & Walden, 2002, P2666, k=0)
r=(2*b+1)/g;
A=sqrt(pi*g*2^r*exp(-gammaln(r)));

%% Morse Wavelet

for i = 1:num_freq

    a = Wbg/(frex(i)*2*pi);
    fa  = a*freq;        
    wa = (2*pi*fa)';
    %wa = [wa,zeros(1,npnts/2-1)];

    % Normalisation factor for unit energy at each scale
    scale_fac=sqrt(srate*a);


    % Generate normalised, scaled wavelet (Olhede & Walden, 2002, P2666, (10), k=0)
    awt = scale_fac*sqrt(2)*A*wa.^b.*(exp(-wa.^g));
    awt = [awt,zeros(1,npnts/2-1)];

    %Plotting

    figure(15), hold on
    plot(hz,awt);
    title("Freq Domain Morse Wavelet: β = " + b + " & γ = " + g)
    xlabel("Hz")
end
