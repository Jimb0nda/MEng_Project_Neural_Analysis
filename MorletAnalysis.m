function ti = waveletAnalysis(dat,trigger, duration, channel, srate)
%% Setup Parameters 

time = -1:1/srate:1;  % best practice to have 0 in the center of the wavelet
st = trigger;
dur = duration;

% Frequency parameters
min_freq = 2; % Hz
max_freq = 30; % Hz
num_freq = 40; % count
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% select channel no
chan_no=channel;

%% variable number of wavelet cycles

% different wavelet widths (number of wavelet cycles)
range_cycles = [ 4 10 ];

% other wavelet parameters
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);

%% Initialise loops


for trial_no = 1:length(st)

    % Indexing for extension phase
    trig_ind=st(trial_no):st(trial_no)+dur(trial_no)-1;

    data = double(squeeze(dat(trig_ind,chan_no)));
    dataR = reshape(data,1,[]);
    
    % Time Axis Setup
    x = (0:length(data)-1)/srate;
    
    % specify baseline periods for dB-conversion
    baseline_window = [ -500  200 ];


    % convert baseline time into indices
    baseidx = dsearchn(x',baseline_window');

    
    % Initialise the time frequency matrix 
    %(number of frequencies by the length of the time vector)
    tf = zeros(num_freq, length(data));
    
    % Initialise the ITPC matrix 
    %(number of frequencies by the length of the time vector)
    %itpc = zeros(num_freq, length(data));
    
    % Create gaussian envelope
    %fwhm = .5;            % width of the Gaussian in seconds
    %gaus = exp(-4*log(2)*time.^2 / fwhm^2);

    % N's for convolution
    nData = length(dataR);
    nKern = length(time);
    nConv = nData + nKern - 1;
    halfK = floor(nKern/2);

    % FFT for data
    dataX = fft(dataR, nConv);

    for fi=1:num_freq

        % create wavelet and get its FFT
        s = nCycles(fi)/(2*pi*frex(fi));
        cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2)); % Morlet Wavelet

        kernel = fft(cmw, nConv);
        
        % max-value normalize the spectrum of the wavelet
        kernel = kernel ./ max(kernel); 
        
        % Convolve
        as = ifft(dataX.*kernel);
        as = as(halfK+1:end-halfK);
        as = reshape(as, length(data),[]); 

        aspow = abs(as).^2;
        
        %Compute time frequency power
        tf(fi,:) = mean(aspow,2);
        tf = tf + tf(fi,:);
        
        % compute ITPC
        %itpc(fi,:) = abs(mean(exp(1j*angle(as)),2));
        %itpc = itpc + itpc(fi,:);

    end
end

% Average time frequency power matrix over number of trials
tf = tf*(1/length(st));

% db normalization
tfDB = 10*log10( bsxfun(@rdivide, tf, mean(tf(:,baseidx(1):baseidx(2)),2)) );


%% Plotting 

figure(1), clf
%Time Frequency Power Plot
s1 = subplot(121);
contourf(x,frex,tf,40,'linecolor','none')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("Time Frequency Power Plot")
colormap(s1,jet);


% dB Normalised Plot
s2 = subplot(122);
contourf(x,frex,tfDB,40,'linecolor','none')
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title("dB Normalised TF Power Plot")
colormap(s2,jet);



