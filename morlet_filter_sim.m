function [eeg, emg, itpc_eeg, itpc_emg] = morlet_filter_sim(srate,dataR1, dataR2, num_freq, frex)

% best practice to have 0 in the center of the wavelet
wavelt = -1:1/srate:1;  

% different wavelet widths (number of wavelet cycles)
range_cycles = [ 4 10 ];
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);


% Initialise the time frequency matrix 
%(number of frequencies by the length of the time vector)
eeg = zeros(num_freq, length(dataR1));
emg = zeros(num_freq, length(dataR1));
itpc_eeg = zeros(num_freq, length(dataR1));
itpc_emg = zeros(num_freq, length(dataR1));

% % N's for convolution
% nData = length(dataR1);
% nKern = length(wavelt);
% nConv = nData + nKern - 1;
% half_wave = floor(nKern/2);

% FFT for eeg and emg data
dataX_eeg = fft(dataR1);
dataX_emg = fft(dataR2);


for fi=1:num_freq

    % create wavelet and get its FFT
    s = nCycles(fi)/(2*pi*frex(fi));
    % Morlet Wavelet
    cmw = exp(2*1i*pi*frex(fi).*wavelt) .* exp(-wavelt.^2./(2*s^2)); 

    kernel = fft(cmw);
    % max-value normalize the spectrum of the wavelet
    kernel = kernel ./ max(kernel); 
    kernel = [kernel,zeros(1,999)];
    

    % Convolve EEG
    eeg_as = ifft(dataX_eeg.*kernel);
    eeg(fi,:) = eeg(fi,:) + eeg_as;
    % Convolve EMG
    emg_as = ifft(dataX_emg.*kernel);
    emg(fi,:) = emg(fi,:) + emg_as;

    % compute ITPC
    itpc_eeg(fi,:) =  exp(1i*angle(eeg_as));
    itpc_emg(fi,:) =  exp(1i*angle(emg_as));

end

end

