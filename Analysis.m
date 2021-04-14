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
eeg_chan=2;  % Ext EEG
emg_chan=3;  % Ext EMG

%% Setup Parameters
srate = 1000;           % in Hz

eeg_tf = 0;
emg_tf = 0;
coherence = 0;
itpc = 0;
t = 0;

min_freq = 4; % Hz
max_freq = 36; % Hz
num_freq = 35; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

%% Loop through trials


for trial_no = 1:length(st1)

    % Indexing for extension phase
    trig_ind=st1(trial_no):st1(trial_no)+2999;

    % Setting up data vectors from dat file
    eeg_data = double(squeeze(dat(trig_ind,eeg_chan)));
    dataR_eeg = reshape(eeg_data,1,[]);
    emg_data = double(squeeze(dat(trig_ind,emg_chan)));
    dataR_emg = abs(reshape(emg_data,1,[]));
    %dataR_emg = reshape(emg_data,1,[]);

    % Morlet Analysis
    %[eeg, emg, itpc] = morlet_filter(srate, dataR_eeg, dataR_emg, num_freq, frex);
    
    % Morse Analysis
    [eeg, emg,t1, itpc_eeg, itpc_emg] = morse_filter(srate,dataR_eeg, dataR_emg, num_freq, frex);
    
    % Time Frequency Cross Spectrum Equations
    eeg_tf = eeg_tf + abs(eeg.*eeg); 
    emg_tf = emg_tf + abs(emg.*emg);
    coherence = coherence + (eeg.*conj(emg));
    itpc = itpc + exp(1j*angle(eeg.*conj(emg)));

    
end

% Average time frequency power matrix and ITPC over number of trials
eeg_tf = eeg_tf/length(st1);
emg_tf = emg_tf/length(st1);
coherence = coherence/length(st1);
itpc = itpc/length(st1);

%% Test

 % Construct output spectral matrix f.
f(:,:,1)=log10(eeg_tf);
f(:,:,2)=log10(emg_tf);
f(:,:,3)=abs(coherence) .* abs(coherence) ./ (eeg_tf.*emg_tf);
f(:,:,4)=abs(itpc);

%% Plotting 

coi = coi();

% Time Axis Setup
timeAxis = (0:length(eeg_data)-1)/srate;

figure(1), clf

% EEG Time Frequency Power Plot
subplot(221);
contourf(timeAxis,frex,f(:,:,1),40,'linecolor','none')
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("EEG Time Frequency Power Plot, channel: " + eeg_chan)

% EMG Time Frequency Power Plot
subplot(222);
contourf(timeAxis,frex,f(:,:,2),40,'linecolor','none')
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("EMG Time Frequency Power Plot, channel: " + emg_chan)

% Coherence Plot
subplot(2,2,3);
hold on
contourf(timeAxis,frex,f(:,:,3),40,'linecolor','none')
colorbar
plot(timeAxis, coi(:,2), 'r')
axis([0 3 4 36]);
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("Coherence Plot for EEG & EMG channels: " + eeg_chan + " & " + emg_chan)
hold off

% ITPC Plot
subplot(2,2,4);
hold on
contourf(timeAxis,frex,f(:,:,4),40,'linecolor','none')
colorbar
plot(timeAxis, coi(:,2), 'r')
axis([0 3 4 36]);
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("ITPC Plot for EEG & EMG channels: " + eeg_chan + " & " + emg_chan)
hold off

% % Line plot for coherence and ITPC plot matching
% figure(2), clf
% hold on
% plot(timeAxis,f(:,:,3),'r')
% plot(timeAxis,f(:,:,4),'k')
% legend({'Coherence';'ITPC'})
% hold off