%% Testing signal analysis techniques
clear

% Load in signal data
load mv_20.mat

%Columns
% 1 EEG 1 (2 cm from midline)
% 2 EEG 2 (4 cm from midline)
% 3 EMG 1 (EDC muscle, wrist extensor)
% 4 EMG 2 (FCR muscle, wrist flexor)
% 5 Acceleration signal

% select channel no
chan_no=2;  % Ext EEG

%% Setup Parameters
srate = 1000;           % in Hz
wavelt = -2:1/srate:2;  % best practice to have 0 in the center of the wavelet

% Frequency parameters
min_freq = 8; % Hz
max_freq = 30; % Hz
num_freq = 40; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% different wavelet widths (number of wavelet cycles)
range_cycles = [ 4 10 ];
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);

%% Initialise loops


for trial_no = 1:length(st1)

    % Indexing for extension phase
    trig_ind=st1(trial_no):st1(trial_no)+2999;

    data = double(squeeze(dat(trig_ind,chan_no)));
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
        
        %Compute time frequency power
        tf(fi,:) = tf(fi,:) + abs(as).^2;
        
        % compute ITPC
        itpc(fi,:) = itpc(fi,:) + abs(mean(exp(1i*angle(as)),2));

    end
end

% Average time frequency power matrix and ITPC over number of trials
tf = tf*(1/length(st1));

%% Plotting 

% Time Axis Setup
timeAxis = (0:length(data)-1)/srate;

figure(1), clf

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
colormap(s2,hot);
