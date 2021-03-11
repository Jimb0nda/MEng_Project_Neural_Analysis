function [eeg, emg, itpc] = morse_filter(srate,dataR_eeg, dataR_emg)

%% Setup Parameters

% best practice to have 0 in the center of the wavelet
wavelt = -2:1/srate:2;

%Wavelet Parameters
a = 0.501;
gamma = 3;
beta = 27;

% Frequency parameters
min_freq = 0;
num_freq = 40; % count

% Initialise the time frequency matrix 
%(number of frequencies by the length of the time vector)
eeg = zeros(num_freq, length(dataR_eeg));
emg = zeros(num_freq, length(dataR_eeg));
itpc = zeros(num_freq, length(dataR_eeg));

% N's for convolution
nData = length(dataR_eeg);
nKern = length(wavelt);
nConv = nData + nKern - 1;
half_wave = floor(nKern/2);

% FFT for eeg and emg data
dataX_eeg = fft(dataR_eeg, nConv);
dataX_emg = fft(dataR_emg, nConv);

for i = 1:beta
    
    Wbg = (i/gamma)^(1/gamma); % maximum value at the peak frequency

    % Frequency parameters
    max_freq = Wbg;
    frex = linspace(min_freq,max_freq,num_freq);
    
    awt =  a .* ((2*pi.*frex).^i) .* (exp(-(2*pi.*frex).^gamma));
    
    kernel = fft(awt, nConv);
    % max-value normalize the spectrum of the wavelet
    kernel = kernel ./ max(kernel); 
    
     % Convolve EEG
    eeg_as = ifft(dataX_eeg.*kernel);
    eeg_as = eeg_as(half_wave+1:end-half_wave);
    eeg(i,:) = eeg(i,:) + eeg_as;
    % Convolve EMG
    emg_as = ifft(dataX_emg.*kernel);
    emg_as = emg_as(half_wave+1:end-half_wave);
    emg(i,:) = emg(i,:) + emg_as;

    % compute ITPC
    itpc(i,:) =  abs(mean(exp(1i*angle(eeg.*conj(emg)))));
    
end

end

