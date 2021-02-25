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
eeg_chan=1;  % Ext EEG
emg_chan=3;  % Ext EMG

%% Setup Parameters
srate = 1000;           % in Hz
wavelt = -2:1/srate:2;  % best practice to have 0 in the center of the wavelet

% Frequency parameters
min_freq = 10; % Hz
max_freq = 35; % Hz
num_freq = 40; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% different wavelet widths (number of wavelet cycles)
range_cycles = [ 4 10 ];
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);

%% Initialise loops


for trial_no = 1:length(st1)

    % Indexing for extension phase
    trig_ind=st1(trial_no):st1(trial_no)+2999;

    % Setting up data vectors from dat file
    eeg_data = double(squeeze(dat(trig_ind,eeg_chan)));
    dataR_eeg = reshape(eeg_data,1,[]);
    emg_data = double(squeeze(dat(trig_ind,emg_chan)));
    dataR_emg = abs(reshape(emg_data,1,[]));

    % Initialise the time frequency matrix 
    %(number of frequencies by the length of the time vector)
    eeg_tf = zeros(num_freq, length(eeg_data));
    emg_tf = zeros(num_freq, length(emg_data));
    coherence = zeros(num_freq, length(emg_data));
    
    % Initialise the ITPC matrix 
    %(number of frequencies by the length of the time vector)
    %itpc = zeros(num_freq, length(eeg_data));

    % N's for convolution
    nData = length(dataR_eeg);
    nKern = length(wavelt);
    nConv = nData + nKern - 1;
    half_wave = floor(nKern/2);

    % FFT for eeg and emg data
    dataX_eeg = fft(dataR_eeg, nConv);
    dataX_emg = fft(dataR_emg, nConv);


    for fi=1:num_freq

        % create wavelet and get its FFT
        s = nCycles(fi)/(2*pi*frex(fi));
        cmw = exp(2*1i*pi*frex(fi).*wavelt) .* exp(-wavelt.^2./(2*s^2)); % Morlet Wavelet

        kernel = fft(cmw, nConv);
        % max-value normalize the spectrum of the wavelet
        kernel = kernel ./ max(kernel); 
        
        % Convolve EEG
        eeg_as = ifft(dataX_eeg.*kernel);
        eeg_as = eeg_as(half_wave+1:end-half_wave);
        % Convolve EMG
        emg_as = ifft(dataX_emg.*kernel);
        emg_as = emg_as(half_wave+1:end-half_wave);
        
        % Time Frequency Cross Spectrum Equations
        Sxy = eeg_as.*conj(emg_as);  %11
        Sx = abs(eeg_as).^2;         %12
        Sy = abs(emg_as).^2;         %13
        
        %Compute time frequency power plots
        eeg_tf(fi,:) = eeg_tf(fi,:) + Sx; 
        emg_tf(fi,:) = emg_tf(fi,:) + Sy;
        
        % Time Frequency Coherence
        coherence(fi,:) = coherence(fi,:) + (abs(Sxy).^2 ./ Sx.*Sy);
        
        % compute ITPC
        %itpc(fi,:) = itpc(fi,:) + abs(mean(exp(1i*angle(Sxy)),2));

    end
end

% Average time frequency power matrix and ITPC over number of trials
eeg_tf = eeg_tf*(1/length(st1));
emg_tf = emg_tf*(1/length(st1));
coherence = coherence*(1/length(st1));

%% Plotting 

% Time Axis Setup
timeAxis = (0:length(eeg_data)-1)/srate;

figure(1), clf

%Time Frequency Power Plot
subplot(221);
contourf(timeAxis,frex,eeg_tf,40,'linecolor','none')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("EEG Time Frequency Power Plot, channel: " + eeg_chan)

%Time Frequency Power Plot
subplot(222);
contourf(timeAxis,frex,emg_tf,40,'linecolor','none')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("EMG Time Frequency Power Plot, channel: " + emg_chan)


% Coherence Plot
subplot(2,2,[3,4]);
contourf(timeAxis,frex,coherence,40,'linecolor','none')
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("Coherence Plot for EEG & EMG channels: " + eeg_chan + " & " + emg_chan + " , During holding phase st1")

% Coherence Plot
%subplot(2,2,4);
%contourf(timeAxis,frex,itpc,40,'linecolor','none')
%xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("Coherence Plot for EEG & EMG channels: " + eeg_chan + " & " + emg_chan + " , During holding phase st1")
