% Optical simulation
%
% Example:  OOK reception of a DPSK signal using a narrow filter
%
% Try to change what happens by:
%       1) changing the duty cycle 
%       2) changing the bandwidth of the optical filter
%       3) changing the extinction ratio

clear all;
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64; % number of symbols
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

x.oftype = 'gauss';   % optical filter type
x.obw    = 0.6;       % optical filter bandwidth (narrow...)
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = inf;       % electrical filter bandwidth (i.e. no filter)
                      % all bandwidths are normalized to the bitrate
x.slopez = 0;         % post-fiber cumulated slope [ps/nm^2]
x.lambda = lam;       % post-fiber central wavelength [nm]
x.dpost  = 0;         % post-fiber cumulated dispersion [ps/nm]

x.rec    = 'ook';     % receiver type
x.plot   = 'ploteye'; % type of plot
x.color  = 'r-';      % color of plot
% x.print  = {'Rxcur','current'};  % type of print to file. By enabling this
                                   % option the current will be printed to a file called 'Rxcur'

% Karhunen-Loeve parameters

x.ber    = 1e-9;       % target bit error rate at which the
                       % OSNR is evaluated.
x.eta    = 4;          % frequency expansion fatcor (see ber_kl.m)
x.mu     = 3.5;        % time expansion fatcor (see ber_kl.m)
x.osnr   = 12+(-10:5); % signal-to-noise ratios [dB]
x.poln   = 2;          % noise polarizations        
x.saddle = 'yes';      % evaluate the BER by saddlepoint 
                       % approximation (see ber_kl.m)

%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

% Preparing the optical field:
Pavg = 1;
E    = lasersource( Pavg, lam, spac);

for ii = 1:Nch
    % Realization of the electrical signals:
    pat(:,ii)    = pattern('debruijn',ii);
    pat_rx(:,ii) = pat_decoder(pat(:,ii),'nf-dpsk');

    elec(:,ii) = electricsource(pat(:,ii),'nf-dpsk',symbrate,'cosroll',duty,roll);
    % Modulating the field;
    Eopt(:,ii) = mz_modulator(E(:,ii), elec(:,ii));
end

create_field('unique',Eopt,[],struct('power','average'));

figure(1); 
plotfield('x',1,'papa','r-')    % plot the Tx field (normalized)
                                % use 'pa--' if you don't want
                                % normalization                             
%%%%%%%% Receiver

figure(2)   % activate figure(3) for the eye diagram
grid on;
hold on;

[pb,osnr]=ber_kl(1,x,pat_rx(:,1));

fprintf('\n\n=========== Results ===========\n\n');
fprintf('OSNR @ BER = %.e : %.4f [dB/0.1 nm]\n',x.ber,osnr);

figure(3)
semilogy(x.osnr,pb);
grid on;
xlabel('OSNR   [dB/0.1nm]')
ylabel('BER')

