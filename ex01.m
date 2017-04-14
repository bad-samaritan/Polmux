% Optical simulation
%
% Example: propagation of a single channel OOK signal within a 100 km of 
%            a fully compensated fiber. 
%
% Try to change what happens by:
%       1) changing the duty cycle (RZ -> duty = 50%)
%       2) pulse roll-off in ELECTRICSOURCE
%       3) Cumulated nonlinear phase

clear all
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64;  % number of symbols
Nt    = 64; % points x symbol
Nch   = 1;   % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi      = 0.4*pi; % peak cumulated nonlinear phase
exratio  = 13;     % extinction ratio [dB]
lam      = 1550;   % central wavelength [nm]
spac     = 0.4;    % channel spacing [nm]
symbrate = 10;     % baud rate [Gbaud]
duty     = 1;      % duty cycle
roll     = 0.2;    % pulse roll-off

%%%%%%%%%%%%%%%% Link parameters

Din   = 0; % in-line dispersion x span [ps/nm]
Nspan = 1; % number of spans

%%%%  Fiber 1 (Tx)
tx.length  = 1e5;     % length [m]
tx.alphadB = 0.2;     % attenuation [dB/km]
tx.aeff    = 80;      % effective area [um^2]
tx.n2      = 2.7e-20; % nonlinear index
tx.lambda  = 1550;    % wavelength [nm] @ dispersion 
tx.disp    = 17;      % dispersion [ps/nm/km] @ wavelength
tx.slope   = 0;       % slope [ps/nm^2/km] @ wavelength
tx.dphimax = 3E-3;    % maximum nonlinear phase rotation per step
tx.dzmax   = 2E4;     % maximum SSFM step 

%%%%  Fiber 2 (compensating fiber)
comp.alphadB = 0.6;     % attenuation [dB/km]
comp.aeff    = 20;      % effective area [um^2]
comp.n2      = 2.7e-20; % nonlinear index
comp.lambda  = 1550;    % wavelength [nm] @ dispersion 
comp.disp    = -100;    % dispersion [ps/nm/km] @ wavelength
comp.slope   = 0;       % slope [ps/nm^2/km] @ wavelength
comp.dphimax = 3E-3;    % maximum nonlinear phase rotation per step
comp.dzmax   = 2E4;     % maximum SSFM step 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Conversions

gam         = 2*pi*tx.n2./(lam.*tx.aeff)*1e21;                         % gamma index [1/W/m]
comp.length = (Din - tx.disp*tx.length/1e3)/comp.disp*1e3;             % comp. fiber length [m]
Gerbio      = tx.alphadB*tx.length/1e3 + comp.alphadB*comp.length/1e3; % ampli gain [dB]

%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch,mfilename);      % reset all global variables
                                        % and write output to ex01.

pat=pattern('debruijn',1);    % generate the bit-pattern
elec=electricsource(pat,'ook',symbrate,'cosroll',duty,roll);

Ppeak = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan); % get the power that
                                    % gives a peak SPM rotation of phi [rad]
Ein    = lasersource(Ppeak, lam, spac);

% Now generate the electric field using a Mach-Zehnder modulator
Eopt = mz_modulator(Ein, elec,struct('exratio',exratio));

% Associate Eopt to the global GSTATE.FIELDX. Note: create_field and
% reset_all must always be called in each simulation.
create_field('unique',Eopt);

figure(1)
plotfield('x',1,'na--','r-')    % plot the Tx field (normalized)
                                % use 'pa--' if you don't want
                                % normalization

%%%%%%%%%% Optical link
fiber(tx,'g-sx') % Tx fiber
fiber(comp,'g-sx')

% Note the flag: 'g-sx', means with GVD, SPM and XPM (in this case it
% corresponds to 'g-s-', since we are single channel. We are thus
% neglecting PMD

ampliflat(Gerbio,'gain')
 
plotfield('x',1,'na--','b-')    % plot the Rx field (normalized)

% The simulation is done. You can check a summary of the system into the
% file simul_out placed in the directory initialized by reset_all.m
