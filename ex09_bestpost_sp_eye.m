% Optical OOK simulation
%
%  Example: Compare the best post compensating dispersion obtained
%           through the functions BEST_EYE and BEST_SP
%

clear all
clc
%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 32; % number of symbols
Nt    = 32; % points x symbol
Nch   = 1;  % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi      = 0.8*pi; % average cumulated nonlinear phase
exratio    = 13;   % extinction ratio [dB]
lam      = 1550;   % central wavelength [nm]
spac     = 0.4;    % channel spacing [nm]
symbrate = 10;     % baudrate [Gbaud]
duty     = 1;      % duty cycle
roll     = 0.2;    % pulse roll-off

%%%%%%%%%%%%%%%% Link parameters

Din   = 0;  % in-line dispersion x span [ps/nm]
Nspan = 20; % number of spans

%%%%  Fiber 1 (Tx)
tx.length  = 1e5;     % length [m]
tx.alphadB = 0.2;     % attenuation [dB/km]
tx.aeff    = 63;      % effective area [um^2]
tx.n2      = 2.7e-20; % nonlinear index
tx.lambda  = 1550;    % wavelength [nm] @ dispersion 
tx.disp    = 8;       % dispersion [ps/nm/km] @ wavelength
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

%%%%%%%%%%%%%%%%% Receiver

x.oftype = 'gauss';   % optical filter type
x.obw    = 1.8;       % optical filter bandwidth 
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = 0.65;      % electrical filter bandwidth
% all bandwidths are normalized to the bitrate

%%% Post-compensating fiber parameters
x.slopez = 0;         % post-fiber cumulated slope [ps/nm^2]
x.lambda = lam;       % post-fiber central wavelength [nm]
x.dpost = [300 1300]; % post-fiber cumulated dispersion [ps/nm]

x.rec = 'ook';        % receiver type

% Karhunen-Lo�ve parameters

x.ber    = 1e-9;       % target bit error rate at which the optimal
                       % post dispersion is searched.
x.eta    = 1.4;        % frequency expansion fatcor (see BER_KL)
x.mu     = 3.5;        % time expansion fatcor (see BER_KL)
x.osnr   = 12+(-1:10); % signal-to-noise ratios [dB]. 
x.poln   = 2;          % noise polarizations        
x.saddle = 'y';        % evaluate the BER by saddlepoint 
                       % approximation (see BER_KL)

% Note that the OSNR in b2b is 12.1 dB, hence we search the OSNR penalty
% within the range -1:10 around 12 dB. By this way we are able to measure
% penalty up to 9.9 dB

% To gain feeling with x.eta and x.mu try to change them.
% Try also to test the accuracy of the saddlepoint approximation by
% removing the line x.saddle = 'y'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Conversions

gam         = 2*pi*tx.n2./(lam.*tx.aeff)*1e21;                         % gamma index [1/W/m]
comp.length = (Din - tx.disp*tx.length/1e3)/comp.disp*1e3;             % comp. fiber length [m]
Gerbio      = tx.alphadB*tx.length/1e3 + comp.alphadB*comp.length/1e3; % ampli gain [dB]

%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

Pavg = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan); % get the average power
E    = lasersource(Pavg, lam, spac);

for ii=1:Nch
    pat(:,ii)=pattern('debruijn',ii);
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii) = mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt,[],struct('power','average'));

%%%%%%%%%% Optical link
for k=1:Nspan
	fiber(tx,'g-sx')
	fiber(comp,'g-sx')
	
	ampliflat(Gerbio,'gain')
end

%%%%%%%%%% Receiver


%%%%%% now search for the optimal post compensation

[be,bpost] = best_sp(1,x,pat);

fprintf('\n ===== Results from best_sp ==========\n\n');
fprintf('\nBest OSNR penalty = %f [dB]\n',be);
fprintf('Best post cumulated dispersion = %f [ps/nm]\n\n',bpost);

% Let us now check the best post dispersion that maximizes the eye opening

[be,bpost] = best_eye(1,x,pat);

fprintf('\n ===== Results from best_eye ==========\n\n');
fprintf('\nBest eye closure penalty = %f [dB]\n',be);
fprintf('Best post cumulated dispersion = %f [ps/nm]\n\n',bpost);

% The two answers are different: optimizing the eye opening is not the same
% as optimizing the BER. The first is clearly a better answer, the second
% is faster.

% Please note that neither best_sp or best_eye change GSTATE.FIELDX. They
% operate over a copy of it.
