function [eeg, emg] = morse_filter(srate,dataR_eeg, dataR_emg,num_freq,frex)

%% Setup Parameters

% best practice to have 0 in the center of the wavelet
wavelt = -2:1/srate:2;

b = 9;
g = 3;

Wbg = (b/g)^(1/g); % maximum value at the peak frequency

% Normalisation factor (Olhede & Walden, 2002, P2666, k=0)
r=(2*b+1)/g;
A=sqrt(pi*g*2^r*exp(-gammaln(r)));

% Initialise the time frequency matrix 
%(number of frequencies by the length of the time vector)
eeg = zeros(num_freq, length(dataR_eeg));
emg = zeros(num_freq, length(dataR_eeg));
%itpc = zeros(num_freq, length(dataR_eeg));

% N's for convolution
nData = length(dataR_eeg);
nKern = length(wavelt);
nConv = nData + nKern - 1;
half_wave = floor(nKern/2);

% FFT for eeg and emg data
dataX_eeg = fft(dataR_eeg, nConv);
dataX_emg = fft(dataR_emg, nConv);

for i = 1:num_freq
    
    y  = Wbg*frex;             % shifted frequencies
    wa = 2*pi*y;

    % Normalisation factor for unit energy at each scale
    scale_fac=sqrt(srate*Wbg);

    % Generate normalised, scaled wavelet (Olhede & Walden, 2002, P2666, (10), k=0)
    awt = scale_fac*sqrt(2)*A*(wa.^b).*(exp(-wa.^g));
    % Botched way of adding the zero padding, would fix later
    awt = ifft(awt);
    awt = fft(awt,nConv);
    % max-value normalize the spectrum of the wavelet
    kernel = awt ./ max(awt); 
    
     % Convolve EEG
    eeg_as = ifft(dataX_eeg.*kernel);
    eeg_as = eeg_as(half_wave+1:end-half_wave);
    eeg(i,:) = eeg(i,:) + eeg_as;
    % Convolve EMG
    emg_as = ifft(dataX_emg.*kernel);
    emg_as = emg_as(half_wave+1:end-half_wave);
    emg(i,:) = emg(i,:) + emg_as;

    % compute ITPC
    %itpc(i,:) =  abs(mean(exp(1i*angle(eeg.*conj(emg)))));
    
end

end

