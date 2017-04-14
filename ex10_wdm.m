% Optical simulation
%
% Example: Same as EX03_WDM, but the channels are propagated into 
%          separate fields.
%
% Note the presence of 'sepfields' in CREATE_FIELD. By this way we 
% propagate 5 separate channels and hence GSTATE.FIELDX has size [N,NCH] 
% being N=GSTATE.NSYMB * GSTATE.NT and NCH=5. Such method allows to neglect 
% FWM and to save computational time since there is no need to use 
% GSTATE.NT large enough to surround all channels. GSTATE.NT can be chosen 
% basing on single channel considerations. Another feature is that it is 
% possible to use the flag '---x' in FIBER, i.e. it is possible to 
% suppress SPM and leave XPM on. Compare the results with EX3_WDM.

clear all
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64; % number of symbols
Nt    = 32; % points x symbol (less point than with unique field)
Nch   = 5;  % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi      = 0.2*pi; % average cumulated nonlinear phase
exratio    = 13;     % extinction ratio [dB]
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

for ii = 1:Nch
    pat(:,ii)=pattern('debruijn',ii); % note the different de Bruijn seeds
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii)=mz_modulator(E(:,ii), elec(:,ii),struct('exratio',exratio));
end

create_field('sepfields',Eopt,[],struct('power','average'));

for k=1:Nch
    figure(k)
    plotfield('x',k,'pa--','r-');   
end

%%%%%%%%%% Optical link

fiber(tx,'g-sx') % Tx fiber
fiber(comp,'g-sx')

% Note the flag: 'g-sx', means with GVD, SPM and XPM. We are thus
% neglecting PMD.

ampliflat(Gerbio,'gain')

for k=1:Nch
    figure(k)
    hold on;
    plotfield('x',k,'pa--','b-'); % plotfield('x',k,'pa--','b-',oftype,obw);
end

% You can compare the results with ex03_wdm. Such a comparison can be done
% by adding the optical filter for a temporal extraction in plotfield, thus
% by uncommenting % plotfield('x',k,'pa--','b-',oftype,obw); above.
% The difference between ex03_wdm and this example is FWM, here negligible.
% Try to decrease the fiber dispersion to observe a significant FWM.
