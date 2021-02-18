clear 

%% Setup Parameters
load sampleEEGdata.mat

% wavelet parameters
num_freq = 40;
min_freq =  2;
max_freq = 30;

channel2use = 'pz';

% set range for variable number of wavelet cycles
range_cycles = [ 3 10 ];

% parameters (notice using logarithmically spaced frequencies!)
frex  = logspace(log10(min_freq),log10(max_freq),num_freq);
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);
wavelt  = -2:1/EEG.srate:2;
half_wave = (length(wavelt)-1)/2;

% FFT parameters
nWave = length(wavelt);
nData = EEG.pnts*EEG.trials;
nConv = nWave+nData-1;


% FFT of data (doesn't change on frequency iteration)
dataX = fft( reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);

% initialize output time-frequency data
itpc = zeros(num_freq,EEG.pnts);

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
    as = reshape(as,EEG.pnts,EEG.trials);

    % compute ITPC
    itpc(fi,:) = abs(mean(exp(1i*angle(as)),2));

end


%% Plotting 

figure(1), clf

% ITPC Plot
contourf(EEG.times,frex,itpc,40,'linecolor','none')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("ITPC")
colormap(hot);
