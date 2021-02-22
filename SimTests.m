clear 

%% Setup Parameters
srate = 1000;           % in Hz
wavelt = -2:1/srate:2;  % best practice to have 0 in the center of the wavelet

% Frequency parameters
min_freq = 4; % Hz
max_freq = 30; % Hz
num_freq = 30; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% different wavelet widths (number of wavelet cycles)
range_cycles = [ 4 10 ];
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);

% time vector
pnts   = 3000;
times = (0:pnts-1)/srate;

for trials = 1:40
    %% transient oscillations w/ Gaussian

    % Gaussian and sine parameters
    peaktime = 2; % seconds
    width    = .12;
    sinefreq = 15; % for sine wave

    % create Gaussian taper
    gaus = exp( -(times-peaktime).^2 / (2*width^2) );

    % trial-unique sine wave
    cosw = cos(2*pi*sinefreq*times);

    data = randn(1,3000) + 0.4*(cosw);% .* gaus

    %% Initialise loops

    dataR = reshape(data,1,[]);

    % Initialise the time frequency matrix 
    %(number of frequencies by the length of the time vector)
    tf = zeros(num_freq, length(data));

    % Initialise the ITPC matrix 
    %(number of frequencies by the length of the time vector)
    itpc = zeros(num_freq, length(data));

    % N's for convolution
    nData = length(dataR);
    nKern = length(wavelt);
    nConv = nData + nKern - 1;
    half_wave = floor(nKern/2);

    % FFT for data
    dataX = fft(dataR, nConv);

    for fi=1:num_freq

        % create wavelet and get its FFT
        s = nCycles(fi)/(2*pi*frex(fi));
        cmw = exp(2*1i*pi*frex(fi).*wavelt) .* exp(-wavelt.^2./(2*s^2)); % Morlet Wavelet

        kernel = fft(cmw, nConv);
        % max-value normalize the spectrum of the wavelet
        kernel = kernel ./ max(kernel); 

        % Convolve
        as = ifft(dataX.*kernel);
        as = as(half_wave+1:end-half_wave);
        %as = reshape(as, length(data),[]); 

        %Compute time frequency power
        tf(fi,:) = tf(fi,:) + abs(as).^2;

        % compute ITPC
        itpc(fi,:) = itpc(fi,:) + abs(mean(exp(1i*angle(as)),2));

    end
end

% Average time frequency power matrix and ITPC over number of trials
tf = tf*(1/40);

%% Plotting 

% Time Axis Setup
timeAxis = (0:length(data)-1)/srate;

figure(10), clf

%Time Frequency Power Plot
s1 = subplot(121);
contourf(timeAxis,frex,tf,40,'linecolor','none')
%conofinf('morl',dataR,length(data),data(1:1594),'plot');
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("Time Frequency Power Plot")
colormap(s1,jet);

% ITPC Plot
s2 = subplot(122);
contourf(timeAxis,frex,itpc,40,'linecolor','none')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("ITPC")
colormap(s2,jet);
