clear 

%% Setup Parameters
srate = 1000;           % in Hz
wavelt = -2:1/srate:2;  % best practice to have 0 in the center of the wavelet

% movie time parameter!
refresh_speed = .6; % seconds

% Frequency parameters
min_freq = 4; % Hz
max_freq = 30; % Hz
num_freq = 30; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% different wavelet widths (number of wavelet cycles)
range_cycles = [ 4 10 ];
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);

%% Wavelet

for fi=1:num_freq
    
        % create wavelet and get its FFT
        s = nCycles(fi)/(2*pi*frex(fi));
        cmw = exp(2*1i*pi*frex(fi).*wavelt) .* exp(-wavelt.^2./(2*s^2)); % Morlet Wavelet

        figure(50), clf
        hold on
        plot(wavelt,real(cmw),'b');
        plot(wavelt,imag(cmw),'r--');
        plot(wavelt,exp(-wavelt.^2./(2*s^2)),'k','linew',3);
        legend({'Real';'Imaginary';'Gaussian'})

        pause(refresh_speed)

end

