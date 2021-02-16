% Script mv20_plot1.m, 4/2/21, DMH
% Extracting and plotting trials from EEG-EMG data MV_20


load mv_20

%Columns
% 1 EEG 1 (2 cm from midline)
% 2 EEG 2 (4 cm from midline)
% 3 EMG 1 (EDC muscle, wrist extensor)
% 4 EMG 2 (FCR muscle, wrist flexor)
% 5 Acceleration signal


% select trials
trial_no=1;

% select channel no
chan_no=3;  % Ext EMG

% Indexing for extension phase
ext_ind=st0(trial_no):st0(trial_no)+dur0(trial_no)-1;

% Indexing for position holding phase (in extension)
pos_ind=st1(trial_no):st1(trial_no)+dur1(trial_no)-1;

% Extract data
dat_ext=dat(ext_ind,chan_no);
dat_pos=dat(pos_ind,chan_no);

figure
subplot(2,1,1)
plot(dat_ext,'k')
title(['Extending. Chan: ',num2str(chan_no),', Trial: ',num2str(trial_no)])
subplot(2,1,2)
plot(dat_pos,'k')
title(['Extended. Chan: ',num2str(chan_no),', Trial: ',num2str(trial_no)])
