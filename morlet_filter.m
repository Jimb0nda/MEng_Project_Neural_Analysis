function [eeg, emg] = morlet_filter(dataR_eeg, dataR_emg, wavelt, num_freq, frex)

% Initialise the time frequency matrix 
%(number of frequencies by the length of the time vector)
eeg = zeros(num_freq, length(eeg_data));
emg = zeros(num_freq, length(emg_data));

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
    eeg(fi,:) = eeg(fi,:) + eeg_as;
    % Convolve EMG
    emg_as = ifft(dataX_emg.*kernel);
    emg_as = emg_as(half_wave+1:end-half_wave);
    emg(fi,:) = emg(fi,:) + emg_as;


end

