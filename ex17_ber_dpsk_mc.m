% Optical DPSK simulation
%
% Example: Evaluate the bit error rate (BER) of a DPSK signal by Monte
%          Carlo Simulation and check it against the Karhunen-Loeve 
%          algorithm.
%


clear all
% clc
%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 32; % number of symbols
Nt    = 32; % points x symbol
Nch   = 1;  % number of channels


reset_all(Nsymb,Nt,Nch);

global CONSTANTS;            % CONSTANTS is a global structure variable.
CLIGHT  = CONSTANTS.CLIGHT;  % speed of light [m/s]
HPLANCK = CONSTANTS.HPLANCK; % Planck's constant [J*s]

%%%%%%%%%%%%%%%%  Pulse parameters

exratio    = 13;   % extinction ratio [dB]
lam      = 1550; % central wavelength [nm]
spac     = 0.4;  % channel spacing [nm]
symbrate = 10;   % baudrate [Gbaud]
duty     = 1;    % duty cycle
roll     = 0.2;  % pulse roll-off

%%%%%%%%%%%%%%%% Monte Carlo parameters

samp.logic = [0;1]; % binary symbols
samp.thr   = 0;     % optimal threshold.

mc.stop = [0.1 68]; % stop criterion: relative error 0.1 with  
                    % Gaussian confidence 68%.


%%%%%%%%%%%%%%%% Receiver
% Note the sampling time fixed to the center of the bit
x.ts     = 0;         % sampling time
x.oftype = 'gauss';   % optical filter type
x.obw    = 2;         % optical filter bandwidth 
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = 0.65;      % i.e., without electrical filter
x.rec    = 'dpsk';    % receiver type

%%%%%%%%%%%%%%%% Karhunen-Loeve parameters (see ber_kl.m)

x.osnr   = 1:4; % signal-to-noise ratios [dB/0.1nm]
x.eta    = 5;   % frequency expansion factor 
x.mu     = 3.5; % time expansion factor 
x.poln   = 2;   % noise polarizations        
x.saddle = 'n'; % saddlepoint approximation? 

%%%%%%%%%%%%%%%% Amplifier parameters

Pavg   = 1;   % average power [mW]
Gerbio = 10;  % gain [dB]
nampli = 1;   % number of amplifiers
osnrbw = 0.1; % bandwidth for OSNR measurement [nm]

hvdl = -30-10*log10(HPLANCK*CLIGHT/lam*CLIGHT*osnrbw/lam^2*1e18); % conv. factor
nsp  = 10*log10(Pavg) + hvdl - 10*log10(10^(Gerbio/10)-1) - 3 - ...
    10*log10(nampli) - x.osnr; % [dB]
F    = nsp + 3; % noise figure [dB]
                % the one-sided ASE PSD of a single amplifier is 2*nsp*hvdl*(Gerbio-1) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('state',1);   % set random state

% First, create the pattern. We plan to transmit a fixed PRBS
% pattern, hence there is no need to insert the call to
% bit_pattern inside the Monte Carlo cycle

for ii = 1:Nch
    pat(:,ii) = pattern('debruijn',ii);
    pat_rx(:,ii) = pat_decoder(pat(:,ii),'dpsk'); % pattern decoding
end

fprintf('\nResults:\n\n');

for knf = 1:length(F)

    ampli.f   = F(knf);

    % Now we are ready for the Monte Carlo simulation
    cond = true;
    while cond % iterate as soon as the stopping condition is reached
        reset_all(Nsymb,Nt,Nch);

        E    = lasersource(Pavg, lam, spac);

        for ii=1:Nch
            elec(:,ii)=electricsource(pat(:,ii),'bpsk',symbrate,'cosroll',duty,roll);
            Eopt = mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
        end

        create_field('unique',Eopt,[],struct('power','average'));

        ampliflat(-Gerbio,'gain');      % noiseless attenuator!
        ampliflat(Gerbio,'gain',ampli); % noisy amplifier

        [eo,ts,y]=eval_eye(1,x,pat_rx(:,1)); % y are the sampled bits
        pat_hat = samp2pat(x,samp,y);
        [cond,avgber,nruns,stdber]=ber_estimate(pat_hat,pat_rx(:,1),mc);
    end % end while cond

    mcber(knf) = avgber;  % Monte Carlo BER
    mcstd(knf) = stdber;  % Monte Carlo error
    mcnruns(knf) = nruns; % Monte Carlo runs

    fprintf('run %2d: BER = %5.2e  std(BER) = %5.2e  # runs=%11d \n',...
        knf,mcber(knf),mcstd(knf),mcnruns(knf));

end

%%%%%%%%%% Now evaluate the same BER using the KL algorithm

x.b2b = 'b2b';
pb = ber_kl(1,x,pat_rx);

figure(1)
semilogy(x.osnr,pb,'b',x.osnr,mcber,'ro');
grid on; hold on;
legend('KL result','Monte carlo BER')
xlabel('OSNR   [dB/0.1 nm]')
ylabel('BER')
