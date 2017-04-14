% Optical simulation
%
% Example: Propagation of a WDM-5 channels OOK signal within a 100 km of 
%          a fully compensated fiber. 
%
% Note: see ex10_wdm to see how to work with separate fields, i.e.
%       without FWM.
%

clear all
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64;  % number of symbols
Nt    = 128; % points x symbol
Nch   = 5;   % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi      = 0.2*pi; % AVERAGE cumulated nonlinear phase
exratio    = 13;   % extinction ratio [dB]
lam      = 1550;   % central wavelength [nm]
spac     = 0.4;    % channel spacing [nm]
symbrate = 10;     % baudrate [Gbaud]
duty     = 1;      % duty cycle
roll     = 0.2;    % pulse roll-off

% Note: the spectral efficiency is 0.2

%%%%%%%%%%%%%%%% Link parameters

Din   = 0; % in-line dispersion x span [ps/nm]
Nspan = 1; % number of spans

%%%%  Fiber 1 (Tx)
tx.length  = 1e5;     % length [m]
tx.alphadB = 0.2;     % attenuation [dB/km]
tx.aeff    = 80;      % effective area [um^2]
tx.n2      = 2.7e-20; % nonlinear index
tx.lambda  = 1550;    % wavelength [nm] @ dispersion 
tx.disp    = 17;       % dispersion [ps/nm/km] @ wavelength
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

%%%%%%%%%%%%%%%% Receiver parameters

oftype = 'butt6'; % optical filter type
obw    = 2.5;     % optical filter bandwidth 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Conversions

gam         = 2*pi*tx.n2./(lam.*tx.aeff)*1e21;                         % gamma index [1/W/m]
comp.length = (Din - tx.disp*tx.length/1e3)/comp.disp*1e3;             % comp. fiber length [m]
Gerbio      = tx.alphadB*tx.length/1e3 + comp.alphadB*comp.length/1e3; % ampli gain [dB]

%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

Pavg = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan); % get the average power
E    = lasersource(Pavg, lam, spac);
% Why Pavg is average since lasersource generates a constant signal? 
% At this stage average or peak is meaningless. It will be create_field that
% will create fields with an average power equal to Pavg. See next.

for ii = 1:Nch
    pat(:,ii)=pattern('debruijn',ii); % note the different de Bruijn seeds
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii) = mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt,[],struct('power','average'));
% here the set: the fields have average power Pavg [mW], i.e. the one
% measured by a power meter placed before the first optical fiber of the
% link.

for k=1:Nch
    figure(k)
    hold on;
    plotfield('x',k,'pa--','b-',oftype,obw); % Note the temporary  
                                             % extraction of channel k.
end

%%%%%%%%%% Optical link
fiber(tx,'g-sx') % Tx fiber
fiber(comp,'g-sx')

% Note the flag: 'g-sx', means with GVD, SPM and XPM (in this case it
% corresponds to 'g-s-', since we are using a unique field. We are thus
% neglecting PMD

ampliflat(Gerbio,'gain')

for k=1:Nch
    figure(k)
    hold on;
    plotfield('x',k,'pa--','r-',oftype,obw); % Note the temporary  
                                             % extraction of channel k.
end

figure
% here we plot the complete spectrum (note the four wave mixing (FWM)
% peaks outside the information bandwidth)
plotfield('x',1,'--p-');
% now we extract channel 4 ONCE AND DEFINITIVELY.
optfilter(4,oftype,obw)
figure
plotfield('x',1,'papa','r');
% Note that after the extraction in the global variable GSTATE.FIELDX we
% have the equivalent lowpass field of channel 4, hence the spectrum is
% still centered around the frequency zero.
