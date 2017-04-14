% DSP for PolMUX Coherent QPSK signal
% Made from Optilux elements (-> dsp4cohdec.m)
function [RxSamples4D] = DspPdmCohQpsk(RxSamples4D, dspParams, chNum)
	global GSTATE;
	global CONSTANTS;  % CONSTANTS is a global structure variable.
	CLIGHT = CONSTANTS.CLIGHT;  % speed of light [m/s]

	% Digital Clock Recovery
	% This operation is useless for the moment. We are already sampling in the
	% middle of the period. Now that we have compensated for chromatic
	% dispersion we need only the central samples.
	if ~dspParams.workatbaudrate
		RxSamples4D = RxSamples4D(1:2:end,:);
	end

	% Non linear phase rotation
	if dspParams.applynlr
		RxSamples4D    = NLRotation(RxSamples4D, dspParams.nlralpha);
	end

	% Normalization
	peak = 4*sqrt(GSTATE.POWER(chNum)); % 2* -> see receiver_cohmix
	RxSamples4D = RxSamples4D/peak;

	% Polar rotation: None, Fixed or Adaptative
	if dspParams.applypol && (size(RxSamples4D,2)==2)
		switch lower(dspParams.polmethod)
			case 'singlepol'
				r = RxSamples4D(:,1) ./ RxSamples4D(:, 2);
				rKikuchi = mean( r );
				RxSamples4D = rotpolar(RxSamples4D, rKikuchi );
			case 'cma'
				RxSamples4D = cmapolardemux( RxSamples4D, dspParams.cmaparams );
			case 'easi'
				RxSamples4D = easipolardemux( RxSamples4D, dspParams.easiparams );
			case 'combo'
				RxSamples4D = easipolardemux( RxSamples4D, dspParams.easiparams );
				RxSamples4D = cmapolardemux( RxSamples4D, dspParams.cmaparams );
			otherwise
				error('Unknown Polar Rotation method.');
		end
	end

	% Carrier Phase Estimation
	M = 2.^dspParams.modorder;
	navg = dspParams.freqavg;
	if navg
		% Estimating frequency in bell labs style [citation neede]:
		omega = cumsum( vitvit( RxSamples4D .* conj( fastshift(RxSamples4D,1) )...
			, M, M, navg , false) );
		% Cleaning omega to match the circularity:
		closestallowedendpoints = omega(1,:) + round((omega(end,:)-omega(1,:))/2/pi)*2*pi;
		correctionratio = closestallowedendpoints ./ omega(end,:);
		omega=((omega-ones(length(omega),1)*omega(1,:)).*(ones(length(omega),1)...
			*correctionratio))+ones(length(omega),1)*omega(1,:);
		% Demodulating signals:
		sigdemod = RxSamples4D .* fastexp( -omega );
		% Estimating phase using Viterbi and Viterbi method:
		navg	= dspParams.phasavg;
		P		= dspParams.poworder;
		theta	= vitvit( sigdemod, P, M, navg, true );
		if dspParams.modorder > 1
			CarrierPhaseOffSet = +pi/4;
		else
			CarrierPhaseOffSet = 0;
		end
		Carrier = fastexp( -omega - theta + CarrierPhaseOffSet );
	else
		navg    = dspParams.phasavg;
		P       = dspParams.poworder;
		theta   = vitvit( RxSamples4D, P, M, navg, true );
		if dspParams.modorder > 1
			CarrierPhaseOffSet = +pi/4;
		else
			CarrierPhaseOffSet = 0;
		end
		Carrier = fastexp(- theta + CarrierPhaseOffSet );
	end
	RxSamples4D = RxSamples4D .*  Carrier ;

	% You can return received and processed signal (as amplitudes and phases)
	%signalPh  = angle(RxSamples4D);
	%signalAmp = abs(RxSamples4D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = NLRotation( s, alpha)
	Phases     = angle( s );
	Amplitudes = abs( s );
	Asquare    = sum( Amplitudes .* Amplitudes , 2 );
	Delta_P    = Asquare - mean( Asquare );
	Phases     = Phases + alpha*Delta_P*ones(1,size(s,2));
	y          = Amplitudes .* fastexp( Phases );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = vitvit( s, P, M, k, applyunwrap)
	L = length(s);
	% Phase -> Phase * M
	% Amplitude -> Amplitude ^ P
	if P == M
		s = s .^ P;
	else
		s = abs(s).^P .* fastexp( angle( s.^ M ) );
	end
	if k>0
		N = 2*k + 1;
		if N<L
			Smoothing_Filter = fft( ones(N, 1) / N , L ) * ones(1,size(s,2));
			s = ifft( fft(s,L) .* Smoothing_Filter );
		else
			slong = repmat( s, ceil(N / L), 1 );
			Smoothing_Filter = fft( ones(N, 1) / N , ceil(N / L).*L ) * ones(1,size(s,2));
			slong = ifft( fft(slong) .* Smoothing_Filter );
			s = slong(1:L,:);
		end
	end
	if applyunwrap
		theta = unwrap( angle(s) ) / M;
	else
		theta = angle( s ) / M;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = rotpolar( x , r )
	if abs(r) < 0.5
		m = abs(r);
		delta = angle(r);
		alpha = m^2 / (m^2 + 1);
	else
		m = abs( 1 / r);
		delta = -angle( 1 / r );
		alpha = 1 / ( m^2 + 1 );
	end
	M = [sqrt(alpha)*exp(-1i*delta) -sqrt(1-alpha)*exp(-1i*delta);...
		 sqrt(1-alpha)				sqrt(alpha)];
	y = x * M;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=cmapolardemux(x, params)
	R       = params.R;			% Radius for convergence
	mu      = params.mu;		% Convergence parameter
	taps    = params.taps;		% length of adaptive filters
	halftaps= floor( taps / 2);	% half of (taps-1)
	polmux  = (params.txpolars == 2); % 1 o 2 polarization transmitted
	hzero = zeros( taps , 2, 2);
	if isfield(params,'mat')
		M = params.mat;
	else
		if ~polmux
			r = mean( x(:,1) ./ x(:,2) );
			M = rotpolar(1, r).';
		else
			phi = params.phizero;
			M = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
		end
	end
	hzero( halftaps+1, :, :) = M;  % Initializing central taps
	if halftaps
		extendedx = [ x(end-halftaps+1:end,:) ; x ; x(1:halftaps,:) ];
	else
		extendedx = x;
	end
	h1 = squeeze( hzero(:, 1, :) );
	h2 = squeeze( hzero(:, 2, :) );
	if taps==1
		h1 = h1.';
		h2 = h2.';
	end
	convergence = false;
	c=1;
	L=length(x);
	repetitions = 50*ceil(1./(L.*mu));
	while ~convergence && (c<(repetitions))
		h1_old = h1*1;
		h2_old = h2*1;
		[y h1_new h2_new] = cmaadaptivefilter(extendedx, h1, h2, taps, mu, R, 1);
		% Here I'm fixing to 1 the number of samples per bit. This can be 2 or
		% more, but syill is not immplemented in dsp4cohdec. In case this
		% changes, the 1 must be replaced with the correct value.
		if any(any(h1_new)) || any(any(h2_new))
			h1=h1_new;
			h2=h2_new;
		end
		if max(max(abs([h1_old-h1 h2_old-h2]))) < 5e-5 % max(max(abs([(h1_old-h1)./(h1+h1_old)/2 (h2_old-h2)./(h1+h1_old)/2])))
			convergence = true;
		end
		c=c+1;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=easipolardemux(x, params)
	mu      = params.mu;	% Convergence parameter
	taps    = 1;             
	halftaps= floor( taps / 2);	% half of (taps-1)
	polmux  = (params.txpolars == 2);	% 1 o 2 polarization transmitted
	hzero = zeros( taps , 2, 2);
	if isfield(params,'mat')
		M = params.mat;
	else
		if ~polmux
			r = mean( x(:,1) ./ x(:,2) );
			M = rotpolar(1, r).';
		else
			phi = params.phizero;
			M = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
		end
	end
	hzero( halftaps+1, :, :) = M;  % Initializing central taps
	if halftaps
		extendedx = [ x(end-halftaps+1:end,:) ; x ; x(1:halftaps,:) ];
	else
		extendedx = x;
	end
	h1 = squeeze( hzero(:, 1, :) );
	h2 = squeeze( hzero(:, 2, :) );
	if taps==1
		h1 = h1.';
		h2 = h2.';
	end
	convergence = false;
	c=1;
	L=length(x);
	repetitions = 20*ceil(1./(L.*mu));
	while ~convergence && (c<(repetitions))
		h1_old = h1*1;
		h2_old = h2*1;
		[y h1_new h2_new] = easiadaptivefilter(extendedx, h1, h2, taps, mu, 1);
		% Here I'm fixing to 1 the number of samples per bit. This can be 2 or
		% more, but syill is not immplemented in dsp4cohdec. In case this
		% changes, the 1 must be replaced with the correct value.
		if any(any(h1_new)) || any(any(h2_new))
			h1=h1_new;
			h2=h2_new;
		end
		if max(max(abs([h1_old-h1 h2_old-h2]))) < 5e-5 % max(max(abs([(h1_old-h1)./(h1+h1_old)/2 (h2_old-h2)./(h1+h1_old)/2])))
			convergence = true;
		end
		c=c+1;
	end
end
