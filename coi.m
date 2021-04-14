function [coi] = coi()
% function [coi] = cwt_coigen(dat_pts,rate,offset_sec,f0,scale_fac)
%
% Function to calculate COI for CWT transform
%
% Inputs
%   dat_pts     No of data points in segment
%   rate        Sampling rate (samples/sec)
%   offset_sec  Any offset from zero in time axis in seconds, added to time axis
%   f0          CWT f0
%   scale_fac   Scaling factor for COI (morlet: sqrt(2), morse: domain_tprime, ...)
%
%
% Output
%   coi         2 column vector with COI: Time in seconds & COI frequency.

%% COI

b = 9;
g = 3;

Wbg = (b/g)^(1/g); % maximum value at the peak frequency


dat_pts = 3002;
rate = 1000;
offset_sec = 0;
f0 = Wbg/(2*pi);
scale_fac = sqrt(2)*(sqrt(b*g)/Wbg);


dat_pts_2=fix((dat_pts-2)/2);
coi_ind=[(1:dat_pts_2) [(1:mod(dat_pts,2))*(dat_pts_2+1)] (dat_pts_2:-1:1)]';
coi_freq=rate*scale_fac*f0./coi_ind;
coi=[(1/rate)*(1:dat_pts-2)'+offset_sec coi_freq];


end
