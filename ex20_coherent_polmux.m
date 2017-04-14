% Coherent PDM-QPSK optical simulation
%
% Example: Propagation of a single channel Polarization Multiplexed QPSK
%          transmission, received  by means of a coherent receiver 
%          featuring digital signal processing
%

clear all
clc

rand('state',1);
randn('state',2);

%%%%%%%%%%%%%%%% Monte Carlo parameters

samp.logic = [0;1]; % binary symbols
samp.thr   = 0;     % optimal threshold.

mc.stop = [0.1 68]; % stop criterion: relative error 0.1 with  
                    % Gaussian confidence 95%.
% mc.nmin = 100;      % minimum number of errors

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 2^6;  % number of symbols
Nt    = 64;   % points x symbol
Nch   = 1;    % number of channels

reset_all(Nsymb,Nt,Nch);

global CONSTANTS;            % CONSTANTS is a global structure variable.
CLIGHT  = CONSTANTS.CLIGHT;  % speed of light [m/s]
HPLANCK = CONSTANTS.HPLANCK; % Planck's constant [J*s]

%%%%%%%%%%%%%%%%%  Pulse parameters
Pavg     = 2;      % total transmitted power [mW]
exratio  = inf;    % extinction ratio [dB]
lam      = 1550;   % central wavelength [nm]
spac     = 0.4;    % channel spacing [nm]
symbrate = 10;     % symbol rate [Gbaud]
duty     = 1.0;    % duty cycle
roll     = 0.2;    % pulse roll-off

%%%%%%%%%%%%%%%%% Receiver
x.rec    = 'coherent'; % receiver type
x.ts     = 0;           
x.oftype = 'gauss';    % optical filter type
x.obw    = 1.9;        % optical filter bandwidth
x.eftype = 'bessel5';  % electrical filter type
x.ebw    = 0.65;       % electrical filter bandwidth
x.osnr   = 7:1:11;     % OSNR [dB/0.1 nm]
                       % all bandwidths are normalized to the symbol-rate
x.delay  = 'theory';

% Coherent Detector Parameters
x.lopower    = 0;    % dBm
% x.lodetuning = 100E6; % 100 MHz (will be automatically rounded to the
                      % nearest admitted frequency)
    
%%%%%%%%%%%%%%%%% Digital Signal Processing Parameters

%   --- General Parameters ---
dspParameters.sps            = Nt;    % samples per received symbol
dspParameters.workatbaudrate = false; % work with one (true) or two (false) 
                                      % samples/symb.
%   --- Analog-to-Digital Converter ---
dspParameters.applyadc     = false;      % apply ADC (true/false)?
dspParameters.adcbits      = 5;          % bits of resolution of the ADC
dspParameters.samplingrate = 2*symbrate; % [GSamples/symbol] sampling rate
%   --- Dispersion Compensating Filter ---
dspParameters.applydcf   = false; % apply DCF (true/false) ?
dspParameters.dispersion = 4000;  % [ps/nm] dispersion to compensate for
dspParameters.ndispsym   = 16;    % length of the filter [symbols]
dspParameters.lambda     = 1550;
%   --- Nonlinear Mitigation ---
dspParameters.applynlr = false;
%   --- Polarization Demultiplexing
dspParameters.applypol  = false;   % apply poldemux (true/false)?
dspParameters.polmethod = 'cma'; % poldemux algorithm

dspParameters.cmaparams.R         = [1 1];  % radius of CMA 
dspParameters.cmaparams.mu        = 1/6000; % convergence parameter of CMA
dspParameters.cmaparams.taps      = 7;      % CMA filter taps
dspParameters.cmaparams.txpolars  = 2;      % transmitted polarizations
dspParameters.cmaparams.phizero   = 0;      % rotation angle between 
                                            % tx and rx field 

dspParameters.easiparams.mu       = 1/6000; % convergence parameter of EASI
dspParameters.easiparams.txpolars = 2;      % transmitted polarizations 
dspParameters.easiparams.phizero  = 0;      % rotation angle between 
                                            % tx and rx field

%   --- Intermediate Frequency Estimation ---
dspParameters.modorder = 2;   % modulation order (BPSK =1, QPSK = 2, etc.)
dspParameters.freqavg  = 500; % phase estimation smoothing parameter
dspParameters.phasavg  = 3;   % phase estimation smoothing parameter (2^N-1)
dspParameters.poworder = 2;   % modulus power in phase estimation. 

%%%%%%%%%%%%%%%%% 

Nmid   = ceil(Nch/2);  % central channel
%%%%%%%%%%%%%%%%% Amplifier parameters
Gerbio = 1; % dummy gain
nampli = 1;                 % number of amplifiers
osnrbw = 0.1;               % bandwidth for OSNR measurement [nm]

hvdl = -30-10*log10(HPLANCK*CLIGHT/lam*CLIGHT*osnrbw/lam^2*1e18); % conv. factor
nsp  = 10*log10(Pavg) + hvdl - 10*log10(10^(Gerbio/10)-1) - 3 - ...
     10*log10(nampli) - x.osnr; % [dB]
F    = nsp + 3; % noise figure [dB]
                % the one-sided ASE PSD of a single amplifier is
                % 2*nsp*hvdl*(Gerbio-1)
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transmission


for ii=1:Nch
    [patx(:,ii) patmatx{ii}]        = pattern('debruijn',1,struct('alphabet',4));
    [paty(:,ii) patmaty{ii}]        = pattern('debruijn',2,struct('alphabet',4));
    [patx_rx(:,ii)  patmatx_rx{ii}] = pat_decoder(patx,'dqpsk');
    [paty_rx(:,ii)  patmaty_rx{ii}] = pat_decoder(paty,'dqpsk');
end

for knf=1:length(F) % cycle over the OSNRs
    fprintf('OSNR = %g \n',x.osnr(knf));

    ampli.f = F(knf); % get amplifier noise figure
    
    cond = true;
    while cond % iterate as soon as the stopping condition is reached
        reset_all(Nsymb,Nt,Nch);

        E = lasersource(Pavg, lam, spac); 
        
        for ii = 1:Nch
            elecx_i(:,ii) = electricsource(patmatx{ii}(:,1),'qpsk',symbrate,'cosroll',duty,roll);
            elecx_q(:,ii) = electricsource(patmatx{ii}(:,2),'qpsk',symbrate,'cosroll',duty,roll);
            elecy_i(:,ii) = electricsource(patmaty{ii}(:,1),'qpsk',symbrate,'cosroll',duty,roll);
            elecy_q(:,ii) = electricsource(patmaty{ii}(:,2),'qpsk',symbrate,'cosroll',duty,roll);
            Eoptx(:,ii)   = qi_modulator(E(:,ii), elecx_i(:,ii), elecx_q(:,ii));
            Eopty(:,ii)   = qi_modulator(E(:,ii), elecy_i(:,ii), elecy_q(:,ii));
        end                             

        create_field('unique',Eoptx,Eopty,struct('power','average'));
        
        ampliflat(-Gerbio,'gain');     
        ampliflat(Gerbio,'gain',ampli);     % noisy amplifier
        
        [phase amplitude] = dsp4cohdec( Nmid, [patx,paty], x, dspParameters);
%         figure(1) % turn off the noise to see a clear constellation
%         polar(phase,amplitude,'ro') % plot the constellation if you like
        
        patmat_hat = samp2pat(x,samp,phase);

        [pat_hat(:,1) patmat_hat(:,1:2)] = pat_decoder(patmat_hat(:,1:2),...
            'dqpsk',struct('binary','true'));        
        [pat_hat(:,2) patmat_hat(:,3:4)] = pat_decoder(patmat_hat(:,3:4),...
            'dqpsk',struct('binary','true'));                

        patmat_rx   = [patmatx_rx{Nmid} patmaty_rx{Nmid}];
        
        % Sometimes it can happen that the two polarizations are switched,
        % i.e. you receive X as Y and viceversa. In this case we switch
        % also the transmitted patterns.
        if sum(sum(patmat_rx(:,1:2)~=patmat_hat(:,3:4))) < sum(sum(patmat_rx(:,1:2)~=patmat_hat(:,1:2)))
            % if comparing tx X and rx X yields more error than comparing
            % tx X with rx Y switch the patterns.
            temp =patmat_hat(:,1:2);
            patmat_hat(:,1:2)=patmat_hat(:,3:4);
            patmat_hat(:,3:4)=temp;
        end
                
        [cond,avgber,nruns,stdber]=ber_estimate(patmat_hat,patmat_rx,mc);
        if mod(nruns,500)
            fprintf('P{error} = %5.2e  std/P{error} = %.3f  # runs=%11d \n',...
                avgber,stdber/avgber,nruns);
        end
                
    end % end while cond
    fprintf('\n');
    mcber(knf) = avgber;    % Monte Carlo BER
    mcstd(knf) = stdber;    % Monte Carlo error
    mcnruns(knf) = nruns;   % Monte Carlo runs

end % end for knf=1:length(F)

figure(1)
semilogy(x.osnr,mcber,'ro');
grid on;
legend('Monte carlo data')
xlabel('OSNR [dB/0.1 nm]')
ylabel('Error probability')
