%
% Example: Propagate an OOK signal in presence of polarization mode
% 		   dispersion (PMD).
%


clear all
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64;  % number of symbols
Nt    = 64; % points x symbol
Nch   = 1;   % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi      = 0.2*pi; % peak cumulated nonlinear phase [rad]
exratio    = Inf;  % extinction ratio [dB]
lam      = 1550;   % central wavelength [nm]
spac     = 0.8;    % channel spacing [nm]
symbrate = 10;     % baudrate [Gbaud]
duty     = 1;      % duty cycle
roll     = 0.2;    % pulse roll-off

%%%%%%%%%%%%%%%% Link parameters

Din   = 0; % in-line dispersion x span [ps/nm]
Nspan = 1; % number of spans

%%%%  Fiber 1 (Tx)
tx.length  = 1e5;     % length [m]
tx.alphadB = 0.2;     % attenuation [dB/km]
tx.aeff    = 63;      % effective area [um^2]
tx.n2      = 2.7e-20; % nonlinear index
tx.lambda  = 1550;    % wavelength [nm] @ dispersion 
tx.disp    = 0;       % dispersion [ps/nm/km] @ wavelength
tx.slope   = 0;       % slope [ps/nm^2/km] @ wavelength
tx.dphimax = 5E-3;    % maximum nonlinear phase rotation per step
tx.dzmax   = 2E4;     % maximum SSFM step 

% ADDITIONAL FLAGS WILL BE SET NEXT

%%%%%%%%%%%%%%%%% Receiver

x.oftype = 'gauss';   % optical filter type
x.obw    = 1.8;       % optical filter bandwidth 
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = 0.65;      % electrical filter bandwidth
                      % all bandwidths are normalized to the bitrate
x.slopez = 0;         % post-fiber cumulated slope [ps/nm^2]
x.lambda = lam;       % post-fiber central wavelength [nm]
x.dpost  = 1000;      % post-fiber cumulated dispersion [ps/nm]
x.plot = 'ploteye';   % plot the eye
x.rec    = 'ook';    % receiver type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Conversions

gam         = 2*pi*tx.n2./(lam.*tx.aeff)*1e21;  % gamma index [1/W/m]
Gerbio      = tx.alphadB*tx.length/1e3; % ampli gain [dB]
Nmid = ceil(Nch/2); % central channel

%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

Ppeak = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan); % get the peak power
E    = lasersource(Ppeak, lam, spac);

pat = pattern('debruijn',1);
elec = electricsource(pat,'ook',symbrate,'cosroll',duty,roll);
Eopt_x = mz_modulator(E, elec,struct('exratio',exratio));
Eopt_y = zeros(Nsymb*Nt,1); % set to zero the y-component

create_field('unique',Eopt_x,Eopt_y); % create x and y

% Note: in absence of Eopt_y Optilux assumes it equal to zero. In such a
%   case a warning is displayed.

figure(1)
grid on; hold on
plotfield('x',1,'p---','k-') % plot the Tx field
title('Tx power')
%%%%%%%%%% Propagation in a PMF fiber

tx.db0 = 0;         % birefringence at GSTATE.FN=0
tx.theta = pi/4;    % waveplates azimuth
tx.epsilon = pi/4;  % waveplates ellipticity
tx.dgd = 0.5; % [symbols] differential group delay (DGD)
tx.nplates = 20;   % number of waveplates (i.e. the fiber is discretized with
	           % tx.nplates PMF fibers, each of length tx.length/tx.nplates

% Since db0,theta,epsilon are fixed the fiber is a PMF where all waveplates
% have the same birefringence. 

fiber(tx,'gp--')
ampliflat(Gerbio,'gain')

figure(1)
plotfield('xy',1,'p---','b-') % plot the Rx field, x and y components
title('Rx power')

% x.dgd=0.5 bits is clearly visible in the received field

% Now check the overall power
figure(2)
grid on; hold on;

plotfield('tot',1,'p---')
title('Total power: PMF case')

% now try to use a random birefringence. We should remove the fixed
% birefringence used in tx:

tx = rmfield(tx,{'db0','theta','epsilon'});

% ok, repeat the measurement
fiber(tx,'gp--')
ampliflat(Gerbio,'gain')

figure(3)
grid on; hold on;
plotfield('tot',1,'p---','g') % note the difference, due to the random 
                              % behavior of the PMD.
title('Total power: random PMD case')                              

%%%%%%%%%% Receiver
figure(4)
pat_rx = pat_decoder(pat,x.rec); % pattern decoding
[eo,ts]=eval_eye(Nmid,x,pat_rx);

fprintf('\n\n=========== Results ===========\n\n');
fprintf('Normalized Eye opening = %.4f\n',eo);
fprintf('Best sampling time = %f\n',ts);


