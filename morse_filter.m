function [eeg, emg] = morse_filter(srate,dataR_eeg, dataR_emg,num_freq,frex)

%% Setup Parameters

% Define complete frequency range
freq_all=[0:length(dataR_eeg)/2,-length(dataR_eeg)/2+1:-1]'*(srate/length(dataR_eeg));
% Apply heaviside function, only positive frequencies, first (T/2)+1 pts.
freq=freq_all(find(freq_all>=0));

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

% FFT for eeg and emg data
dataX_eeg = fft(dataR_eeg);
dataX_emg = fft(dataR_emg);

for i = 1:num_freq
    
    a = Wbg/(frex(i)*2*pi);
    y  = a*freq;             % shifted frequencies
    wa = (2*pi*y)';
    %wa = [wa,zeros(1,length(dataR_eeg)/2-1)];

    % Normalisation factor for unit energy at each scale
    scale_fac=sqrt(srate*a);

    % Generate normalised, scaled wavelet (Olhede & Walden, 2002, P2666, (10), k=0)
    awt = scale_fac*sqrt(2)*A*(wa.^b).*(exp(-wa.^g)); 
    awt = [awt,zeros(1,length(dataR_eeg)/2-1)];
    
    
    
     % Convolve EEG
    eeg_as = ifft(dataX_eeg.*awt);
    eeg(i,:) = eeg(i,:) + eeg_as;
    % Convolve EMG
    emg_as = ifft(dataX_emg.*awt);
    emg(i,:) = emg(i,:) + emg_as;

    % compute ITPC
    %itpc(i,:) =  abs(mean(exp(1i*angle(eeg.*conj(emg)))));
    
end

end

