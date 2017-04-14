% Optical simulation
%
% Example: propagation of a single channel OOK signal within a 100 km of 
%          a fully compensated fiber. Based on the GUI fiber for an 
%          interactive plot.
%

clear all
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nbit = 32;        % number of bits
Nt = 32;          % points x bit
Nch = 1;          % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi = 0.2*pi;     % peak cumulated nonlinear phase
exratio = 13;     % extinction ratio [dB]
lam = 1550;       % central wavelength [nm]
spac=0.4;         % channel spacing [nm]
symbrate = 10;    % bit rate [Gb/s]
duty = 1;         % duty cycle
roll = 0.2;       % pulse roll-off

%%%%%%%%%%%%%%%% Link parameters

Din = 0;        % in-line dispersion x span [ps/nm]
Nspan = 1;      % number of spans

%%%%  Fiber 1 (Tx)
Lf = 100000;        % length [m]
alphadB = 0.2;      % attenuation [dB/km]
Aeff = 80;          % effective area [um^2]
n2= 2.7e-20;        % nonlinear index
lambda = 1550;      % wavelength [nm] @ dispersion 
Dc = 17;            % dispersion [ps/nm/km] @ wavelength
Slope = 0;          % slope [ps/nm^2/km] @ wavelength

%%%%  Fiber 2 (compensating fiber)
alphadB2 = 0.6;
Aeff2 = 20;
n22= 2.7e-20;
lambda2 = 1550;
Dc2 = -100;
Slope2 = 0;

%%%% GUI options
                            % 1: flag on. 0: flag off.
infoax.flag1d = [1 1 0];    % [Power Phase Chirp]: plot power/phase
infoax.flag3d = [0 0 0];    % [Power Phase Chirp]: 3D not activated
infoax.ch = 1;              % channel number
infoax.comp = 1;            % local full compensation. Each step, a virtual 
                            % optical fiber is inserted to recover the
                            % dispersion. Such option allows to clearly
                            % separate the GVD impairment to the nonlinear
                            % one.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Conversions

gam = 2*pi*n2./(lam.*Aeff)*1e21;            % gamma index [1/W/m]
Lf2 = (Din - Dc*Lf/1e3)/Dc2*1e3;            % comp. fiber length [m]
Gerbio = alphadB*Lf/1e3 + alphadB2*Lf2/1e3; % ampli gain [dB]

%%%%%%%%%% Tx side

reset_all(Nbit,Nt,Nch);

Pavg = phi2pow(phi,Lf,alphadB,gam,0,Nspan); % get the average power
E    = lasersource(Pavg, lam, spac);

pat=pattern('debruijn',1);

elec=electricsource(pat,'ook',symbrate,'cosroll',duty,roll);
Eopt = mz_modulator(E,elec,struct('exratio',exratio));

create_field('unique',Eopt);

%%%%%%%%%% Optical link
fibergui(Lf,alphadB,Aeff,n2,lambda,Dc,Slope,20000,0.005,'gsx',infoax) % Tx fiber
fibergui(Lf2,alphadB2,Aeff2,n22,lambda,Dc2,Slope2,30000,0.005,'gsx',infoax)

% Note the flag: 'gsx', means with GVD, SPM and XPM (in this case it
% corresponds to 'gs-', since we are single channel.

ampliflat(Gerbio,'gain')
 
% Operate on the dynamic interactive figure
