% Optical DPSK simulation
%
% Example: Evaluate the bit error rate (BER) of a DPSK signal and check
%          it against the exact expression for ideal system with matched  
%          optical filter.
%


clear all
clc
%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 32; % number of symbols
Nt    = 32; % points x symbol
Nch   = 1;  % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

exratio    = 13;   % extinction ratio [dB]
lam      = 1550; % central wavelength [nm]
spac     = 0.4;  % channel spacing [nm]
symbrate = 10;   % baudrate [Gbaud]
duty     = 1;    % duty cycle
roll     = 0.2;  % pulse roll-off


%%%%%%%%%%%%%%%%% Receiver

x.oftype = 'movavg'; % optical filter type (ideal in time domain)
x.obw    = 1;        % optical filter bandwidth -> moving average
                     % of 1 bit time
x.eftype = 'gauss';  % electrical filter type
x.ebw    = Inf;      % i.e., without electrical filter
                     % all bandwidths are normalized to the bitrate
x.rec    = 'dpsk';   % receiver type

% Karhunen-Lo�ve parameters

x.ber    = 1e-5;      % target bit error rate at which the
                      % OSNR is evaluated.
x.eta    = 15;        % frequency expansion fatcor (see ber_kl.m)
x.mu     = 5.5;       % time expansion factor (see ber_kl.m)
x.osnr   = 3+(-1:10); % signal-to-noise ratios [dB]
x.poln   = 1;         % noise polarizations        
x.saddle = 'n';       % evaluate the BER without saddlepoint 
                      % approximation? (see ber_kl.m)
x.delay  = 'theory';

reset_all(Nsymb,Nt,Nch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The BER of a DPSK using a matched filter is [see Proakis, Digital
% Communications, McGraw Hill, 4th ed]]:
%
% BER=1/2*exp(-Eb/N0)
%
% where Eb/N0 is the signal-to-noise ratio over a bandwidth equal to the
% bitrate and using one noise polarization. Let us check this result with
% the Karhunen Loeve method implemented in BER_KL.
%
% Note that it is    Eb/N0 = OSNR*dv
%
% being dv=clight/lambda^2*B/bitrate, where B=0.1 nm is the reference 
% bandwidth of the OSNR in BER_KL.

global CONSTANTS;  % CONSTANTS is a global structure variable.
CLIGHT = CONSTANTS.CLIGHT; % speed of light [m/s]

EbN0 = 10.^(x.osnr/10)*CLIGHT/symbrate/lam^2*0.1; % Eb/N0 factor 
EbN0dB = 10*log10(EbN0);

%%%%%%%%%% Tx side

Pavg = 1;
E    = lasersource(Pavg, lam, spac);

for ii=1:Nch
    pat(:,ii)=pattern('11'); % ideal condition of the formula
    elec(:,ii)=electricsource(pat(:,ii),'bpsk',symbrate,'cosroll',duty,roll);
    Eopt(:,ii)=mz_modulator(E(:,ii), elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt);

figure(1)
plotfield('x',1,'na--','r-') % plot the Tx field (normalized to Tx peak 
                             % power).


%%%%%%%%%% Receiver
pat_rx = pat_decoder(pat,'dpsk'); % pattern decoding
[pb,osnr_b2b]=ber_kl(1,x,pat_rx);

pbexact = 1/2*exp(-EbN0);

figure(2)
semilogy(EbN0dB,pb,'b',EbN0dB,pbexact,'r-o');
grid on; hold on;
legend('KL result','Exact BER')
xlabel('Eb/N0   [dB]')
ylabel('BER')
