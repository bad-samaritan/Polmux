% My PDM-QPSK simulation

%%% Initialize

%close all;
clear all;
clc;

rand('state',1);
randn('state',2);

if 'Signal parameters'
	sigParams.logic = [0;1];	% Binary symbols
	sigParams.thr   = 0;		% Optimal threshold
	Pavg = 2;			% Total transmitted power [mW]
	chSpac = 0.4;		% Channel spacing [nm]
	symbolRate = 10;	% Symbol rate [Gbaud]
	duty = 1.0;			% Duty cycle
	roll = 0.2;			% Pulse roll-off,  0.2 = almost square
	% Tx Bits number = 4 * Nsymb //QPSK
	Nsymb = 4096;	% Number of symbols
	Nt	  = 64;		% Points per symbol (optical samples)
	Nch	  = 1;		% Number of channels
	lam	  = 1310;	% Central wavelength [nm]
	OSNR  = 16;		% Set OSNR
	reset_all(Nsymb,Nt,Nch);
end

if 'Fiber parameters'
	fiberOn = true;
	fib.length  = 1e3;		% Length [m]
	fib.alphadB = 0.2;		% Attenuation [dB/km]
	fib.aeff    = 80;		% Effective area [um^2]
	fib.n2      = 2.7e-20;	% Nonlinear index
	fib.lambda  = lam;		% Wavelength [nm] @ dispersion 
	% fib.disp = 60;  NOT OK
	% fib.disp = 100; OK
	fib.disp    = 60;		% Dispersion [ps/nm/km] @ wavelength
	fib.slope   = 0;		% Slope [ps/nm^2/km] @ wavelength
	fib.dphimax = 5E-3;		% Maximum nonlinear phase rotation per step
	fib.dzmax   = 2E4;		% Maximum SSFM step
	%%% Fiber PMD
	fib.pmd		= 0.04;		% Average PMD value ps/sqrt(km)
	fib.nplates = 10;		% Number of waveplates
	fib.manakov = 'no';		% Solve Manakov equation istead of CNLSE
	fib.dgd		= fib.pmd*sqrt(fib.length/1e3)*symbolRate*10-3;
							% Effective DGD [symbols]
end

if 'Receiver parameters'
	%%%% Filters settings
	RxParams.rec    = 'coherent';	% Receiver type
	RxParams.ts     = 0;			% Fixed sampling time
	RxParams.oftype = 'gauss';		% Optical filter type
	RxParams.obw    = 1.9;			% Optical filter bandwidth
	RxParams.oord   = 3;			% Optical filter order (from VPI.TxFilter)
	RxParams.eftype = 'bessel5';	% Electrical filter type
	RxParams.ebw    = 0.65;			% Electrical filter bandwidth
	RxParams.eord   = 4;			% Electrical filter order (from VPI.Rx)
	RxParams.delay  = 'theory';		% Delay from GSTATE.Delay
	RxParams.lopower = 0;			% Coherent Detector Parameters (dBm)
	%%%% ADC parameters
	RxParams.sps            = Nt;		% Samples per received symbol
	RxParams.workatbaudrate = false;	% Work with one (true) or two (false) samples/symb.
	RxParams.applyadc		= true;		% Apply ADC (true/false)?
	RxParams.adcbits     	= 5;		% Bits of resolution of the ADC
	RxParams.baudrate		= symbolRate;	% Signal baudrate
	RxParams.samplingrate	= 2*symbolRate;	% [GSamples/symbol] Sampling rate
	%%%% Build-in Dispersion Compensating Filter
	RxParams.applydcf   = false;	% Apply DCF (true/false) ?
	RxParams.dispersion = 17e5;		% [ps/nm] Dispersion to compensate for
	RxParams.ndispsym   = 16;		% Length of the filter [symbols]
	RxParams.lambda     = lam;		% Operating lambda
end

if 'DSP parameters'
	dspParams.workatbaudrate = false; % Work with one (true) or two (false) samples/symb.
	%%%% Nonlinear Mitigation
	dspParams.applynlr = false;
	%%%% Polarization Demultiplexing
	dspParams.applypol  = false;	% apply poldemux (true/false)?
	dspParams.polmethod = 'cma';	% poldemux algorithm
	dspParams.cmaparams.R         = [1 1];  % radius of CMA 
	dspParams.cmaparams.mu        = 1/6000; % convergence parameter of CMA
	dspParams.cmaparams.taps      = 7;      % CMA filter taps
	dspParams.cmaparams.txpolars  = 2;      % transmitted polarizations
	dspParams.cmaparams.phizero   = 0;      % rotation angle between tx and rx field 
	dspParams.easiparams.mu       = 1/6000; % convergence parameter of EASI
	dspParams.easiparams.txpolars = 2;      % transmitted polarizations 
	dspParams.easiparams.phizero  = 0;      % rotation angle between tx and rx field
	%%%% Intermediate Frequency Estimation
	dspParams.modorder = 2;   % modulation order (BPSK =1, QPSK = 2, etc.)
	dspParams.freqavg  = 500; % phase estimation smoothing parameter
	dspParams.phasavg  = 3;   % phase estimation smoothing parameter (2^N-1)
	dspParams.poworder = 2;   % modulus power in phase estimation. 
end

fprintf('######## %s ########\n', char(datetime));
%%% Sending data
reset_all(Nsymb,Nt,Nch);
carrierElec = lasersource(Pavg, lam, chSpac);
for chNum = 1:Nch
	% Generete PRBS(debruijn) and QPSK symbols
	[TxSymbolsX, TxBitsX] = pattern('debruijn',chNum+1,struct('alphabet',4));
	[TxSymbolsY, TxBitsY] = pattern('debruijn',chNum+2,struct('alphabet',4));
	TxBits4D{chNum} = [TxBitsX TxBitsY];
	TxSymbolsXY{chNum} = [TxSymbolsX TxSymbolsY];
	
	% Generate electric field
	ElecIX(:,chNum) = electricsource(TxBits4D{chNum}(:,1),'qpsk',symbolRate,'cosroll',duty,roll);
	ElecQX(:,chNum) = electricsource(TxBits4D{chNum}(:,2),'qpsk',symbolRate,'cosroll',duty,roll);
	ElecIY(:,chNum) = electricsource(TxBits4D{chNum}(:,3),'qpsk',symbolRate,'cosroll',duty,roll);
	ElecQY(:,chNum) = electricsource(TxBits4D{chNum}(:,4),'qpsk',symbolRate,'cosroll',duty,roll);
	OptSignalX(:,chNum) = qi_modulator(carrierElec(:,chNum), ElecIX(:,chNum), ElecQX(:,chNum));
	OptSignalY(:,chNum) = qi_modulator(carrierElec(:,chNum), ElecIY(:,chNum), ElecQY(:,chNum));
end
create_field('sepfields', OptSignalX, OptSignalY, struct('power','average'));

%%% Transmission
if fiberOn
	% On-Off flags: g - GVD, p - PMD, s - SPM, x - XPM 
	fiber(fib,'g---');
end

%%% Receiving data
for chNum = 1:Nch
	% Receive and process
	[RxSamplCompXY{chNum}, worsteyeop] = RxPdmCohQpsk(...
		chNum, [TxSymbolsXY{chNum}(:,1) TxSymbolsXY{chNum}(:,2)], RxParams);
	% Plot unprocessed signal
	figure(100+chNum);
	polar(angle(RxSamplCompXY{chNum}), abs(RxSamplCompXY{chNum}), 'y.');
	hold on;
	polar(angle(RxSamplCompXY{chNum}(1:2:end,:)), abs(RxSamplCompXY{chNum}(1:2:end,:)), 'g.');
	title('Received samples'); %legend({'Non-Decimated','Decimated'});
	hold off;
	
	% Plot processed but uncompensated signal
	[OutSampCompNonCDXY{chNum}] = DspPdmCohQpsk(RxSamplCompXY{chNum}, dspParams, chNum);
	figure(200+chNum);
	polar(angle(OutSampCompNonCDXY{chNum}), abs(OutSampCompNonCDXY{chNum}), 'r.');
	title('Unprocessed samples');
	
	% Process with custom CD compensation algorithms
	if 'Custom CD comp'
		samplingRateIn = RxParams.samplingrate * 1e9;
		lambdaRef = lam * 1e-9;
		span = fib.length;
		D = fib.disp * 1e-6;
		S = fib.slope * 1e-6;
		fftLength = 256;
		L = fftLength/2;
		NTaps = 300;
		[RxSamplCompCD_XY{chNum}(:,1), RxSamplCompCD_XY{chNum}(:,2), ~] = ...
			CDE_OFDE(RxSamplCompXY{chNum}(:,1), RxSamplCompXY{chNum}(:,2),...
			samplingRateIn, lambdaRef, span, D, S, fftLength, L);
		% Plot CD compensated signal
		figure(300+chNum);
		RxSamplPhDecimCD{chNum} = angle(RxSamplCompCD_XY{chNum});
		RxSamplAmpDecimCD{chNum}  = abs(RxSamplCompCD_XY{chNum});
		polar(angle(RxSamplCompCD_XY{chNum}), abs(RxSamplCompCD_XY{chNum}), 'c.');
		hold on;
		polar(angle(RxSamplCompCD_XY{chNum}(1:2:end,:)), abs(RxSamplCompCD_XY{chNum}(1:2:end,:)), 'b.');
		title('CD compensated samples');
		hold off;
	else
		RxSamplCompCD_XY = RxSamplCompXY;
	end
	
	% Process signal with DSP
	[OutSampCompXY{chNum}] = DspPdmCohQpsk(RxSamplCompCD_XY{chNum}, dspParams, chNum);
	RxSignalPhase{chNum} = angle(OutSampCompXY{chNum});
	RxSignalAmpl{chNum}  = abs(OutSampCompXY{chNum});
	% Plot processed signal
	figure(400+chNum);
	polar(RxSignalPhase{chNum}, RxSignalAmpl{chNum}, 'k.');
	title('Processed samples');
	
	% Decode QPSK
	RxBits4D{chNum} = samp2pat(RxParams, sigParams, RxSignalPhase{chNum});
end

% Calc Errors
for chNum = 1:Nch
	rBitsNoX{chNum} = length(TxBits4D{chNum}(:,1)) + length(TxBits4D{chNum}(:, 2));
	mBitsNoX{chNum} = sum(sum(TxBits4D{chNum}(:,1:2) == RxBits4D{chNum}(:,1:2)));
	eBitsNoX{chNum} = rBitsNoX{chNum} - mBitsNoX{chNum};
	rBitsNoY{chNum} = length(TxBits4D{chNum}(:,1)) + length(TxBits4D{chNum}(:, 2));
	mBitsNoY{chNum} = sum(sum(TxBits4D{chNum}(:,1:2) == RxBits4D{chNum}(:,1:2)));
	eBitsNoY{chNum} = rBitsNoY{chNum} - mBitsNoY{chNum};
	fprintf('Ch %d Pol X Match: %d / %d | Errors: %d\n', chNum, mBitsNoX{chNum}, rBitsNoX{chNum}, eBitsNoX{chNum});
	fprintf('Ch %d Pol Y Match: %d / %d | Errors: %d\n', chNum, mBitsNoY{chNum}, rBitsNoY{chNum}, eBitsNoY{chNum});
end



