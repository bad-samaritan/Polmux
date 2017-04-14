% Optical DQPSK simulation
%
% Example: Evaluate the bit error rate (BER) of a DQPSK signal by Monte
%          Carlo Simulation and check it against the Karhunen-Lo�ve 
%          algorithm.
%


clear all
clc
%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64; % number of symbols
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

mc.stop = [0.1 99]; % stop criterion: relative error 0.1 with  
                    % Gaussian confidence 95%.
% mc.nmin = 100;       % minimum number of errors


%%%%%%%%%%%%%%%% Receiver
% Note the sampling time fixed to the center of the bit
x.ts     = 0;         % sampling time
x.oftype = 'gauss';   % optical filter type
x.obw    = 2;         % optical filter bandwidth 
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = 0.65;      % i.e., without electrical filter
x.rec    = 'dqpsk';    % receiver type
x.comp='both';

%%%%%%%%%%%%%%%% Karhunen-Loeve parameters (see ber_kl.m)

x.osnr   = 2:6; % signal-to-noise ratios [dB/0.1nm]
x.eta    = 5;    % frequency expansion factor 
x.mu     = 3.5;  % time expansion factor 
x.poln   = 2;    % noise polarizations        
x.saddle = 'n';  % saddlepoint approximation? 
%%%%%%%%%%%%%%%% Amplifier parameters

Pavg   = 1;   % average power [mW]
Gerbio = 10;  % [dB]
nampli = 1;   % number of amplifiers
osnrbw = 0.1; % bandwidth for OSNR measurement [nm]

hvdl = -30-10*log10(HPLANCK*CLIGHT/lam*CLIGHT*osnrbw/lam^2*1e18); % conv. factor
nsp  = 10*log10(Pavg) + hvdl - 10*log10(10^(Gerbio/10)-1) - 3 - ...
    10*log10(nampli) - x.osnr; % [dB]
F    = nsp + 3; % noise figure [dB]
                % the one-sided ASE PSD of a single amplifier is 2*nsp*hvdl*(Gerbio-1) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('state',2);   % set random state

% First, create the pattern. We plan to transmit a fixed PRBS
% pattern, hence there is no need to insert the call to
% bit_pattern inside the Monte Carlo cycle

for ii = 1:Nch
    [pat(:,ii) patmat{ii}] = pattern('debruijn',ii,struct('alphabet',4));
    [pat_rx(:,ii) patmat_rx{ii}] = pat_decoder(pat(:,ii),'dqpsk');
end


fprintf('\nResults:\n\n');

% Now we are ready for the Monte Carlo simulation
for knf=1:length(F)  % cycle over the OSNRs

    ampli.f = F(knf);

    cond = true;
    while cond % iterate as soon as the stopping condition is reached
        reset_all(Nsymb,Nt,Nch);       
        
        % Preparing the optical field:
        E    = lasersource(Pavg, lam, spac);

        for ii = 1:Nch
            elec1(:,ii) = electricsource(patmat{ii}(:,1),'qpsk',symbrate,'cosroll',1,0.2);
            elec2(:,ii) = electricsource(patmat{ii}(:,2),'qpsk',symbrate,'cosroll',1,0.2);            
            Eopt(:,ii)  = qi_modulator(E(:,ii), elec1, elec2);
        end
        
        create_field('unique',Eopt,[],struct('power','average'));
        
        ampliflat(-Gerbio,'gain');          % noiseless attenuator!
        ampliflat(Gerbio,'gain',ampli);     % noisy amplifier
        
        [eo,ts,y]=eval_eye(1,x,patmat_rx{ii}); % y are the sampled bits
        pat_hat = samp2pat(x,samp,y);
        [cond,avgber,nruns,stdber]=ber_estimate(pat_hat,patmat_rx{1},mc);
        
    end % end while cond

    mcber(knf) = avgber;    % Monte Carlo BER
    mcstd(knf) = stdber;    % Monte Carlo error
    mcnruns(knf) = nruns;   % Monte Carlo runs

    fprintf('OSNR %2g: BER = %5.2e  std(BER) = %5.2e  # runs=%11d \n',...
        x.osnr(knf),mcber(knf),mcstd(knf),mcnruns(knf));
end % end for knf=1:length(F)

%%%%%%%%%% Now evaluate the same BER using the KL algorithm

x.b2b = 'b2b';
x.comp = 'phase';
pb_i = ber_kl(1,x,patmat_rx{1}(:,1));
x.comp = 'quadrature';
pb_q = ber_kl(1,x,patmat_rx{1}(:,2));

figure(1)
semilogy(x.osnr,pb_i,'b',x.osnr,pb_q,'r',x.osnr,mcber,'bo');
grid on;
legend('KL result I','KL result Q','Monte carlo BER')
xlabel('OSNR [dB/0.1 nm]')
ylabel('BER')
