% 4-point PolMUX Coherent QPSK Receiver
% Made from Optilux elements (-> dsp4cohdec.m)
function [RxSamples, worsteyeop] = RxPdmCohQpsk(chNum, symbolPattern, RxParams)
	global GSTATE;
	global CONSTANTS;  % CONSTANTS is a global structure variable.
	CLIGHT = CONSTANTS.CLIGHT;  % speed of light [m/s]
	
	% Photodiodes Order:
	% Channel 1 : Polar 1, In phase
	% Channel 2 : Polar 2, In phase
	% Channel 3 : Polar 1, In quadrature
	% Channel 4 : Polar 2, In quadrature
	COS_POL1 = 1;
	COS_POL2 = 3;
	SIN_POL1 = 2;
	SIN_POL2 = 4;

	switch RxParams.rec
		case 'coherent'        
			[Irxt,RxParams] = receiver_cohmix(chNum,RxParams); % get the current
			% if GSTATE.FIELDY is created along propagation, it doesn't carry
			% information, thus we must ignore it and consider only the X
			% polarization for which the pattern is defined
		otherwise
			error('Flag X.rec must be ''coherent''');
	end
	if isempty(GSTATE.FIELDY) || (size(symbolPattern,2) == 1)
		Irx = Irxt(:,1:2); % use X polarization
		RxParams.isy = 0;
	else
		Irx = Irxt;        % use both polarizations
		RxParams.isy = 1;
	end

	% ADC finite bit number emulator
	if RxParams.applyadc
		M = max( max( abs( Irx ) ) );   % <- Optimal range. There's no loss.
		% Quantization:
		Irx = round( ( Irx + M ) / 2 ./ M * 2^RxParams.adcbits )* 2 .* M / 2^RxParams.adcbits - M;
	end
	[eyeb,best_ts,delay]=mygeteyeinfo(chNum,Irx,RxParams,symbolPattern);
	for npol = 1:size(Irx,2)/2
		Irx(:,(1:2)+2*(npol-1)) = fastshift(Irx(:,(1:2)+2*(npol-1)),round(-delay(npol)*GSTATE.NT));
	end
	% we will sample in the middle of the period

	% The ADC works at a given sampling rate. This aspect cannot be neglected
	% like the finite bit number. So now there is the mandatory decimation:
	if RxParams.workatbaudrate
		DecimationRate = RxParams.sps; % We will have one sample x symbol
	else
		DecimationRate = RxParams.sps / 2; % We will have two samples x symbol
	end
	if size( Irx,2) == 2
		Irxdec = [ decimate( Irx(:,1), DecimationRate,16,'fir' ) ...
				decimate( Irx(:,2), DecimationRate,16,'fir' ) ];
	else
		Irxdec = [ decimate( Irx(:,1), DecimationRate,16,'fir' ) ...
				decimate( Irx(:,2), DecimationRate,16,'fir' ) ...
				decimate( Irx(:,3), DecimationRate,16,'fir' ) ...
				decimate( Irx(:,4), DecimationRate,16,'fir' ) ];
	end
	% Rebuilding Signals from Sequences ____________________________________
	switch size( Irxdec, 2 ) 
		case 2
			RxSamples = complex( Irxdec(:,1), Irxdec(:,2) );
		case 4
			RxSamples = [ complex( Irxdec(:,COS_POL1), Irxdec(:,SIN_POL1) ) ...
						complex( Irxdec(:,COS_POL2), Irxdec(:,SIN_POL2) ) ];
		otherwise
			error('Unable to manage Irx matrix.');
	end

	% Digital Dispersion Compensation______________________________________
	if RxParams.applydcf
		Beta2L           = -RxParams.dispersion * RxParams.lambda^2 / 2 / pi / CLIGHT * 1E-21;
		SamplesPerSymbol = 1 + ~RxParams.workatbaudrate;
		DCF_Points       = length(RxSamples);
		DCF_Length       = RxParams.ndispsym*SamplesPerSymbol;
		DCF_BandW        = SamplesPerSymbol*RxParams.baudrate; % bitrate in Hz
		Hfilt            = DispCompFilter( Beta2L, DCF_BandW, DCF_Points , DCF_Length );
		Hfilt			 = Hfilt * ones(1,size(RxSamples,2));
		RxSamples		 = ifft( fft(RxSamples) .* Hfilt );
	end

	worsteyeop = min(eyeb(abs(eyeb-mod(eyeb,pi/2))<1e-10));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hfilt = DispCompFilter( Beta2L, B , N, FilterLength)
	freq = ( -B/2 : B / N : B/2 * (N-2) / N );
	freq = ifftshift( freq );
	delay = 2*pi*freq/B * ( FilterLength/2  );
	argum = (2*pi*freq).^2 * Beta2L / 2 - delay;
	H = fastexp ( argum );
	b = ifft( H );
	b = b(1:FilterLength + 1 );
	Hfilt = fft(b.', N) .* fastexp(delay.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eyeb,best_ts,delay,xopt]=mygeteyeinfo(ich,Iric,x,pat)
	% Eye evaluations. On input:
	%   ICH: channel number
	%   IRIC: electric current returned by the receiver
	%   X: struct variable (see EVAL_EYE)
	%   PAT: symbols pattern
	%
	%   On output:
	%
	%   EYEB: Eye opening
	%   BEST_TS: best normalized sampling time (center of bit is 0)
	%   DELAY: overal normalized (to the bit time) delay
	%   XOPT: Best sampling time in discrete points.
	%

	global GSTATE

	for npol = 1:1+x.isy
		ipat = pat(:,npol);

		ipat_old = ipat;
		ipat(ipat == 0) = -3/4*pi;
		ipat(ipat == 1) = 3/4*pi;
		ipat(ipat == 2) = -1/4*pi;
		ipat(ipat == 3) = 1/4*pi;

		%         Nfft = GSTATE.NSYMB*GSTATE.NT;
		%         refsig = reshape(repmat(ipat,1,GSTATE.NT)',Nfft,1);
		if isfield(x,'delay') && strcmp(x.delay,'theory') % theoretical delay
			if isfield(x,'b2b')
				if strcmp(x.b2b,'b2b'), avgdelay = 0;end   % back-to-back: no delay
			else
				if x.isy
					avgdelay = 0.5*(GSTATE.DELAY(1,ich)+GSTATE.DELAY(2,ich));
				else
					avgdelay = GSTATE.DELAY(1,ich);
				end
			end
			delay = ones(1,1+x.isy)*(avgdelay+evaldelay(x.oftype,x.obw*0.5)+...
				evaldelay(x.eftype,x.ebw)+x.post_delay);
			Iric_t = angle(complex(Iric(:,1+2*(npol-1)),...
				Iric(:,2+2*(npol-1))));
		else
			[delay(npol),wrn,rho,Iric_t] = corrdelay(complex(Iric(:,1+2*(npol-1)),...
				Iric(:,2+2*(npol-1))),ipat,GSTATE.NT,GSTATE.NSYMB,'phase');
		end
		halfbit = GSTATE.NT/2;

		nshift = round(halfbit-delay(npol)*GSTATE.NT); % the first bit is centered at index 1
		Iricmat = reshape(fastshift(Iric_t,nshift),GSTATE.NT,GSTATE.NSYMB*size(Iric_t,2))';  % Note the transpose!

		for nalph = 1:max(ipat_old)+1
			minv(:,nalph+(max(ipat_old)+1)*(npol-1)) = min(Iricmat(ipat_old==nalph-1,:),[],1);
			maxv(:,nalph+(max(ipat_old)+1)*(npol-1)) = max(Iricmat(ipat_old==nalph-1,:),[],1);
		end

		if npol == 1
			Iricmatx = Iricmat;
			Iricx       = Iric_t;
		else
			Iricmaty = Iricmat;
			Iricy       = Iric_t;
		end
	end

	if (~isempty(minv) && ~isempty(maxv))

		eyeop(:,1) = minv(:,3)-maxv(:,1);    % eye opening
		eyeop(:,2) = minv(:,2)-maxv(:,4);
		eyeop(:,3) = minv(:,4)-maxv(:,3);
		eyeop(:,4) = minv(:,1)-(maxv(:,2)-2*pi);
		if x.isy
			eyeop(:,5) = minv(:,7)-maxv(:,5);    % eye opening
			eyeop(:,6) = minv(:,6)-maxv(:,8);
			eyeop(:,7) = minv(:,8)-maxv(:,7);
			eyeop(:,8) = minv(:,5)-(maxv(:,6)-2*pi);
		end
		worsteye = min(eyeop,[],2);

		if isfield(x,'ts') % fixed sampling time
			xopt = round((x.ts+0.5)*GSTATE.NT);
			eyeb = eyeop(xopt,:);
			eyeb(eyeb<0) = NaN;
			best_ts = x.ts;
		else
			[eyeb,best_tsn] = max(worsteye);
			eyeb =eyeop(best_tsn,:);

			if GSTATE.NT == 2
				best_ts = best_tsn/GSTATE.NT-0.5;
			else    % Interpolation for the best eye opening
				if best_tsn == 1
					to = 1:3;
				elseif best_tsn == GSTATE.NT
					to = best_tsn-2:best_tsn;
				else
					to = best_tsn-1:best_tsn+1;
				end
				% Parabolic interpolation around the max using three points.
				xopt = 0.5*(to(3)^2*(worsteye(to(1))-worsteye(to(2)))+to(1)^2*(worsteye(to(2))-worsteye(to(3)))+...
					to(2)^2*(worsteye(to(3))-worsteye(to(1))))/(to(3)*(worsteye(to(1))-worsteye(to(2)))+...
					to(1)*(worsteye(to(2))-worsteye(to(3)))+to(2)*(worsteye(to(3))-worsteye(to(1))));
				best_ts = xopt/GSTATE.NT-0.5;
				for nalph = 1:size(eyeop,2)
					eyeb(nalph) = (xopt-to(2))*(xopt-to(3))/((to(1)-to(2))*(to(1)-to(3)))*eyeop(to(1),nalph)+...
						(xopt-to(1))*(xopt-to(3))/((to(2)-to(1))*(to(2)-to(3)))*eyeop(to(2),nalph)+...
						(xopt-to(1))*(xopt-to(2))/((to(3)-to(1))*(to(3)-to(2)))*eyeop(to(3),nalph);
				end
			end
		end
	else     % if (~isempty(min1) && ~isempty(max0))
		eyeb = Inf;     % the eye does not exist
		best_ts = 0;    % any time is good
		xopt = 0.5*GSTATE.NT;
	end

	if isfield(x,'plot')
		if ~isfield(x,'color')
			col = 'b-';
		else
			col = x.color;
		end

		if  strcmpi(x.plot,'ploteye')
			ntime = -0.5:1/(GSTATE.NT-1):0.5;
			if ~isfield(x,'pol')
				if x.isy, 
					x.pol = 'xy';
				else
					x.pol = 'x';
				end
			end
			if isfield(x,'pol')
				if strcmp(x.pol,'y')
					plot(ntime,Iricmaty',col)
				elseif strcmp(x.pol,'xy')
					subplot(1,2,1)
					plot(ntime,Iricmatx',col);
					subplot(1,2,2)
					plot(ntime,Iricmaty',col)
				else
					plot(ntime,Iricmatx',col)
				end
			end
			xlabel('time')
			ylabel('Eye')
		elseif strcmpi(x.plot,'plotcur')
			time=0:1/GSTATE.NT:GSTATE.NSYMB-1/GSTATE.NT;
			if isfield(x,'pol') && strcmp(x.pol,'y')
				plot(time,Iricy,col);
			else
				plot(time,Iricx,col);
			end
			xlabel('time   [a.u.]');
			ylabel('Current   [mA]')
		elseif strcmp(x.plot,'')
			% nothing to do
		else
			error(['the field plot must be ',...
				'''ploteye'' or ''plotcur''']);
		end
		drawnow;
	end
end
