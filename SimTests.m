clear 

%% Setup Parameters
srate = 1000;           % in Hz
wavelt = -2:1/srate:2;  % best practice to have 0 in the center of the wavelet

data = 0;
sim_tf1 = 0;
sim_tf2 = 0;
coherence = 0;
itpc = 0;

% Frequency parameters
min_freq = 1; % Hz
max_freq = 50; % Hz
num_freq = 30; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% different wavelet widths (number of wavelet cycles)
range_cycles = [ 4 10 ];
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);

% time vector
pnts   = 3000;
times = (0:pnts-1)/srate;

% Impulse tests
data = zeros(1,3000);
data(1) = 100; 

% % Gaussian and sine parameters
% peaktime = [0.5 1 1.5 2 2.5]; % seconds
% width    = .12;
% sinefreq = [5 10 15 20 25]; % for sine wave
% 
% for i = 1:length(peaktime)
%    
%     % create Gaussian taper
%     gaus = exp( -(times-peaktime(i)).^2 / (2*width^2) );
%     % trial-unique sine wave
%     cosw = cos(2*pi*sinefreq(i)*times);
%     
%     data = data + 0.4*(cosw).*gaus;
%     
% end


% figure(9), clf
% subplot(211)
% plot(data)
% title("Cosine Pulses at Increasing Frequency")
% 
% data1 = data  +  noise;
% snr = snr(data,noise);
% 
% subplot(212)
% plot(data1)
% title("Cosine Pulses at Increasing Frequency with added Noise, SNR:" + snr + "dB")

%% Loop

for trials = 1:40
    
    noise1 = randn(1,3000);
    noise2 = randn(1,3000);
    data1 = data + noise1;
    data2 = data + noise2;

    dataR1 = reshape(data1,1,[]);
    dataR2 = reshape(data2,1,[]);

    % Morlet Analysis
    %[tf1, tf2, itpc_eeg, itpc_emg] = morlet_filter(srate, dataR1, dataR2, num_freq, frex);
    
    % Morse Analysis
    [tf1, tf2, itpc_eeg, itpc_emg] = morse_filter(srate,dataR1, dataR2, num_freq, frex);
    
    
    % Time Frequency Cross Spectrum Equations
    sim_tf1 = sim_tf1 + abs(tf1.*tf1); 
    sim_tf2 = sim_tf2 + abs(tf2.*tf2);
    coherence = coherence + (tf1.*conj(tf2));
    % compute ITPC
    itpc = itpc + exp(1j*angle(tf1.*conj(tf2)));

    
end

% Average time frequency power matrix and ITPC over number of trials
sim_tf1 = sim_tf1/40;
sim_tf2 = sim_tf2/40;
coherence = coherence/40;
itpc = itpc/40;

 % Construct output spectral matrix f.
f(:,:,1)=log10(sim_tf1);
f(:,:,2)=log10(sim_tf2);
f(:,:,3)=abs(coherence) .* abs(coherence) ./ (sim_tf1.*sim_tf2);
f(:,:,4)=abs(itpc);

%% Plotting 

coi = coi();

% Time Axis Setup
timeAxis = (0:length(data)-1)/srate;

figure(12), clf

%Time Frequency Power Plot
%s1 = subplot(121);
hold on
contourf(timeAxis,frex,f(:,:,3),40,'linecolor','none')
colorbar
title("Impulse at t = 1.5s, Morse Wavelet \beta = 9 & \gamma = 3")
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(jet);
plot(timeAxis, coi(:,2), 'k')
axis([0 3 1 50]);
hold off

%figure(11), clf
% ITPC Plot
% s2 = subplot(122);
% contourf(timeAxis,frex,f(:,:,4),40,'linecolor','none')
% colorbar
% xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("ITPC")
% colormap(jet);
