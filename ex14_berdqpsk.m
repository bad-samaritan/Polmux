% Optical DQPSK simulation
%
% Example: Evaluate the bit error rate (BER) of a DQPSK signal and check
%          it against the exact expression for ideal system with matched
%          optical filter.
%

clear all;
clc
% close all;


%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64;  % number of symbols
Nt    = 64; % points x symbol
Nch   = 1;   % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

exratio    = 13;   % extinction ratio [dB]
lam      = 1550; % central wavelength [nm]
spac     = 0.4;  % channel spacing [nm]
symbrate = 10;   % baudrate [Gbaud]
duty     = 1;    % duty cycle
roll     = 0.2;  % pulse roll-off

%%%%%%%%%%%%%%%%% Receiver

x.oftype = 'gauss';   % optical filter type
x.obw    = 1;         % optical filter bandwidth
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = inf;       % electrical filter bandwidth
                      % all bandwidths are normalized to the symbol-rate
x.slopez = 0;         % post-fiber cumulated slope [ps/nm^2]
x.lambda = lam;       % post-fiber central wavelength [nm]
x.dpost  = 0;         % post-fiber cumulated dispersion [ps/nm]

x.rec    = 'dqpsk';   % receiver type
x.comp   = 'phase';
x.plot   = 'ploteye'; % type of plot
x.color  = 'r-';      % color of plot
% x.print = {'Rxcur','current'};  % type of print to file. By enabling this
% option the current will be printed to a file called 'Rxcur'

% Karhunen-Loeve parameters

x.ber    = 1e-5;        % target bit error rate at which the
                        % OSNR is evaluated.
x.eta    = 10;          % frequency expansion fatcor (see BER_KL)
x.mu     = 5.5;         % time expansion fatcor (see BER_KL)
x.osnr   = 12+(-10:15); % signal-to-noise ratios [dB]
x.poln   = 1;           % noise polarizations        
x.saddle = 'n';         % evaluate the BER by saddlepoint 
                        % approximation? (see BER_KL)
x.delay  = 'theory';    % Use the theoretical delay, Note that in this very 
                        % special case the signal is constant, hence it is
                        % not possible to measure the delay by correlation.
                        
%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

% Preparing the optical field:
Pavg = 1;
E    = lasersource( Pavg, lam, spac);

% Creating the electrical signals:
for ii=1:Nch
    pat{1}(:,ii) = pattern('00');
    pat{2}(:,ii) = pattern('00');
    pat_tx{ii} = [pat{1}(:,ii) pat{2}(:,ii)];
    [temp pat_rx{ii}] = pat_decoder(pat_tx{ii},'dqpsk',struct('binary',true));
    elec{1}(:,ii)=electricsource(pat{1}(:,ii),'qpsk',symbrate,'cosroll',duty,roll);
    elec{2}(:,ii)=electricsource(pat{2}(:,ii),'qpsk',symbrate,'cosroll',duty,roll);
    Eopt(:,ii) = qi_modulator(E(:,ii), elec{1}(:,ii),elec{2}(:,ii));
end

create_field('unique',Eopt);

figure(2);

global CONSTANTS; % CONSTANTS is a global structure variable.
CLIGHT = CONSTANTS.CLIGHT; % speed of light [m/s]

EbN0 = 10.^(x.osnr/10)*CLIGHT/symbrate/lam^2*0.1/2; % Eb/N0 factor 
EbN0dB = 10*log10(EbN0);

a=sqrt(2*EbN0*(1-cos(0.25*pi)));
b=sqrt(2*EbN0*(1+cos(0.25*pi)));

% pb_exact=marcumq(a,b,1)-0.5*exp(-0.5*(a.^2+b.^2)).*besseli(0,a.*b);
% the following line approximates the above commented in order to not use
% the marcumq Matlab function
pb_exact=besseli(0,a.*b)./exp(a.*b).*(exp(-(b-a).^2/2)...
    +a.*sqrt(pi/2).*erfc((b-a)/sqrt(2)))-0.5*exp(-0.5*(a.^2+b.^2)).*besseli(0,a.*b);

[pb_i osnr] = ber_kl(1,x,pat_rx{ii}(:,1)); % in-phase component
x.comp='quadrature';
[pb_q osnr] = ber_kl(1,x,pat_rx{ii}(:,2)); % quadrature component

pb=(pb_i+pb_q)/2;

figure(4)
semilogy(EbN0dB,pb,'b-o');
grid on; hold on;
semilogy(EbN0dB,pb_exact,'r--d');
legend('KL','exact ber');

% let us take look at the constellation (just one symbol in this case)
global GSTATE;

figure(5);
ang=angle(GSTATE.FIELDX);
rho=abs(GSTATE.FIELDX);
polar(ang,rho,'o');
