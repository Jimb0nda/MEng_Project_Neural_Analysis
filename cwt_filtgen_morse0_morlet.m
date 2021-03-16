function [cwt_filt]=cwt_filtgen_morse0_morlet(samp_rate,seg_pts,params,scale_list);

%function [cwt_filt]=cwt_filtgen_morse0_morlet(samp_rate,seg_pts,cwt_type,params,scale_list);
%
% samp_rate  Sampling rate
% seg_pts    Total no of points in each segment, can include zero padding, must be even number
% params     Structure with wavelet parameters
%             params.cwt_type (string): 'morlet' or 'morse'.
%             params.f0:  wavelet f0 (morlet and morse)
%             params.b, params.g: beta and gamma values (morse only)
% scale_list List of CWT to generate

n_scale=length(scale_list);
if (n_scale==0)
	error('No scales to calculate.')
end	

% Define matrix for CWT filters, includes heaviside H(w): no negative frequency values.
cwt_filt=zeros(seg_pts/2+1,n_scale);

% Define complete frequency range
freq_all=[0:seg_pts/2,-seg_pts/2+1:-1]'*(samp_rate/seg_pts);
% Apply heaviside function, only positive frequencies, first (T/2)+1 pts.
freq=freq_all(find(freq_all>=0));

switch params.cwt_type
	case 'morlet'
		f0=params.f0;
		for ind_s=1:n_scale
			% Specific scale
			a=scale_list(ind_s);
			% Normalisation factor for unit energy at each scale
			scale_fac=sqrt(2*a*samp_rate)*(pi^0.25);
			
			% Generate normalised, scaled wavelet
			cwt_filt(:,ind_s)=scale_fac*exp(-0.5*(2*pi*a*freq-2*pi*f0).^2);
		end

	case 'morse'
		b=params.b;  % beta
		g=params.g;  % gamma
		% Normalisation factor (Olhede & Walden, 2002, P2666, k=0)
		r=(2*b+1)/g;
		A=sqrt(pi*g*2^r*exp(-gammaln(r)));
		for ind_s=1:n_scale
			% Specific scale
			a=scale_list(ind_s);
			fa=a*freq;
			wa=2*pi*fa;
			% Normalisation factor for unit energy at each scale
			scale_fac=sqrt(samp_rate*a);
			
			% Generate normalised, scaled wavelet (Olhede & Walden, 2002, P2666, (10), k=0)
			cwt_filt(:,ind_s)=scale_fac*sqrt(2)*A*wa.^b.*exp(-wa.^g);
		end

	otherwise
		error(['Unknown wavelet type: ',params.cwt_type])
end	
