% Coherent QPSK optical simulation
%
% Example: Propagation of a single channel QPSK transmission, received 
%          by means of a coherent receiver featuring digital signal 
%          processing
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

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 2^8;  % number of symbols
Nt    = 64;   % points x symbol
Nch   = 1;    % number of channels
Npol  = 1;    % number of polarizations

reset_all(Nsymb,Nt,Nch);

global CONSTANTS;            % CONSTANTS is a global structure variable.
CLIGHT  = CONSTANTS.CLIGHT;  % speed of light in vacuum [m/s]
HPLANCK = CONSTANTS.HPLANCK; % Planck's constant [J*s]

%%%%%%%%%%%%%%%%%  Pulse parameters
Pavg     = 1;    % transmitted power per channel [mW]
exratio  = inf;  % extinction ratio [dB]
lam      = 1550; % central wavelength [nm]
spac     = 0.4;  % channel spacing [nm]
symbrate = 10;   % symbol rate [Gbaud]
duty     = 1.0;  % duty cycle
roll     = 0.2;  % pulse roll-off

%%%%%%%%%%%%%%%%% Receiver
x.rec = 'coherent';   % receiver type
x.ts     = 0;
x.oftype = 'gauss';   % optical filter type
x.obw    = 1.9;       % optical filter bandwidth
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = 0.65;      % electrical filter bandwidth
                      % all bandwidths are normalized to the symbol-rate
x.osnr   = 3:0.5:7.5; % OSNR [dB/0.1 nm]

x.delay  = 'theory'; 

% Coherent Detector Parameters
x.lopower    = 0;     % dBm
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
dspParameters.applynlr  = false;
%   --- Polarization Demultiplexing
dspParameters.applypol           = false;    % apply polar demux (true/false ) ?
%   --- Intermediate Frequency Estimation ---
dspParameters.modorder = 2;   % modulation order (BPSK =1, QPSK = 2, etc.)
dspParameters.freqavg  = 500; % phase estimation smoothing parameter
dspParameters.phasavg  = 3;   % phase estimation smoothing parameter (2^N-1)
dspParameters.poworder = 2;   % modulus power in phase estimation. See [2]

%%%%%%%%%%%%%%%%% Init

Nmid   = ceil(Nch/2); % central channel

%%%%%%%%%%%%%%%%% Amplifier parameters

Gerbio = 1;         % dummy gain
nampli = 1;         % number of amplifiers
osnrbw = 0.1;       % bandwidth for OSNR measurement [nm]

hvdl = -30-10*log10(HPLANCK*CLIGHT/lam*CLIGHT*osnrbw/lam^2*1e18); % conv. factor
nsp  = 10*log10(Pavg) + hvdl - 10*log10(10^(Gerbio/10)-1) - 3 - ...
     10*log10(nampli) - x.osnr; % [dB]
F    = nsp + 3; % noise figure [dB]
                % the one-sided ASE PSD of a single amplifier is
                % 2*nsp*hvdl*(Gerbio-1)
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transmission

for ii = 1:Nch
    [pat(:,ii) patmat{ii}]         = pattern('debruijn',1,struct('alphabet',4));
    [pat_rx(:,ii)  patmat_rx{ii} ] = pat_decoder(pat(:,ii),'dqpsk');
end

for knf=1:length(F)  % cycle over the OSNRs
    fprintf('OSNR = %g \n',x.osnr(knf));

    ampli.f = F(knf);
    
    cond = true;
    while cond % iterate as soon as the stopping condition is reached
        reset_all(Nsymb,Nt,Nch);

        % options.linewidth = 1E6 / symbrate /1E9;
        E = lasersource(Pavg, lam, spac);
        
        for ii = 1:Nch
            elec_i(:,ii) = electricsource(patmat{ii}(:,1),'qpsk',symbrate,'cosroll',duty,roll);
            elec_q(:,ii) = electricsource(patmat{ii}(:,2),'qpsk',symbrate,'cosroll',duty,roll);

            Eopt(:,ii) = qi_modulator(E(:,ii), elec_i(:,ii), elec_q(:,ii));
        end
        
        create_field('unique',Eopt,[],struct('power','average'));       
             
        ampliflat(-Gerbio,'gain');      % noiseless amplifier
        ampliflat(Gerbio,'gain',ampli); % noisy amplifier
        % evaluate the symbols phase and amplitude (constellation)
        [phase amplitude] = dsp4cohdec( 1, pat(:,1), x, dspParameters);
%         figure(1) % turn off the noise to see a clear constellation
%         polar(phase,amplitude,'ro') % plot the constellation if you like
        pat_hat = samp2pat(x,samp,phase);
        
        [pat_hat patmat_hat] = pat_decoder(pat_hat,'dqpsk',struct('binary','true'));
        
        [cond,avgber,nruns,stdber]=ber_estimate(patmat_hat,patmat_rx{Nmid},mc);
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
