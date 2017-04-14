% Optical DPSK simulation
%
% Example: propagation of a single channel DPSK signal within a 100 km of 
%          a fully compensated fiber followed by a post-compensating 
%          fiber. 
%   
%       1) Evaluate the eye with a fixed post-compensating fiber dispersion    
%       2) Search the best post compensating fiber that minimizes the eye 
%          closure penalty.
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
tx.aeff    = 80;      % effective area [um^2]
tx.n2      = 2.7e-20; % nonlinear index
tx.lambda  = 1550;    % wavelength [nm] @ dispersion 
tx.disp    = 17;       % dispersion [ps/nm/km] @ wavelength
tx.slope   = 0;       % slope [ps/nm^2/km] @ wavelength
tx.dphimax = 5E-3;    % maximum nonlinear phase rotation per step
tx.dzmax   = 2E4;     % maximum SSFM step 

%%%%  Fiber 2 (compensating fiber)
comp.alphadB = 0.6;     % attenuation [dB/km]
comp.aeff    = 20;      % effective area [um^2]
comp.n2      = 2.7e-20; % nonlinear index
comp.lambda  = 1550;    % wavelength [nm] @ dispersion 
comp.disp    = -100;    % dispersion [ps/nm/km] @ wavelength
comp.slope   = 0;       % slope [ps/nm^2/km] @ wavelength
comp.dphimax = 5E-3;    % maximum nonlinear phase rotation per step
comp.dzmax   = 2E4;     % maximum SSFM step 


%%%%%%%%%%%%%%%%% Receiver

x.oftype = 'gauss';   % optical filter type
x.obw    = 1.8;       % optical filter bandwidth 
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = 0.65;      % electrical filter bandwidth
                      % all bandwidths are normalized to the bitrate
x.slopez = 0;         % post-fiber cumulated slope [ps/nm^2]
x.lambda = lam;       % post-fiber central wavelength [nm]
x.dpost  = 1000;      % post-fiber cumulated dispersion [ps/nm]

x.rec    = 'dpsk';    % receiver type
x.plot   = 'ploteye'; % type of plot
x.color  = 'r-';      % color of plot
% x.print = {'Rxcur','current'};  % type of print to file. By enabling this
% option the current will be printed to a file called 'Rxcur'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Conversions

gam         = 2*pi*tx.n2./(lam.*tx.aeff)*1e21;                         % gamma index [1/W/m]
comp.length = (Din - tx.disp*tx.length/1e3)/comp.disp*1e3;             % comp. fiber length [m]
Gerbio      = tx.alphadB*tx.length/1e3 + comp.alphadB*comp.length/1e3; % ampli gain [dB]
Nmid = ceil(Nch/2); % central channel

%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

Ppeak = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan); % get the peak power
E    = lasersource(Ppeak, lam, spac);

for ii=1:Nch
    pat(:,ii)=pattern('debruijn',ii);
    elec(:,ii)=electricsource(pat(:,ii),'bpsk',symbrate,'cosroll',duty,roll);
    Eopt(:,ii)=mz_modulator(E(:,ii), elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt); % associate Eopt to GSTATE.FIELDX

figure(1)
hold on
plotfield('x',1,'pa--','r-','gauss',3) % plot the Tx field

%%%%%%%%%% Optical link
fiber(tx,'g-sx')
fiber(comp,'g-sx')

ampliflat(Gerbio,'gain')

plotfield('x',1,'pa--','b-','gauss',3) % plot the Rx field of ch.1

%%%%%%%%%% Receiver

figure(2)
pat_rx = pat_decoder(pat(:,Nmid),x.rec); % pattern decoding
[eo,ts]=eval_eye(Nmid,x,pat_rx);

fprintf('\n\n=========== Results ===========\n\n');
fprintf('Normalized Eye opening = %.4f\n',eo);
fprintf('Best sampling time = %f\n',ts);
fprintf('Post fiber dispersion =%.2f [ps/nm]\n',x.dpost);


%%%%%% now search for the optimal post compensation
dpostini = x.dpost;
x.dpost = [-1500 1500];  % range for the search

[be,bpost] = best_eye(Nmid,x,pat_rx);
% During the search of the best post you will see the following warning:
% Warning: Plot or print turned off during the search of the optimal post.

fprintf('\nBest Eye closure penalty = %.4f [dB]\n',be);
fprintf('Best post cumulated dispersion = %.2f [ps/nm]\n\n',bpost);
fprintf('\nOther details into the summary file simul_out\n'); 


figure(2)   % Re-plot the best eye, just for comparison
hold on
title(['RED: post = ',num2str(dpostini),' ps/nm,  BLUE: best post = ',...
    num2str(bpost),' ps/nm']);

x.dpost = bpost;  % use the best post
x.color = 'b-';   % change color
[eo2,ts2]=eval_eye(Nmid,x,pat_rx);

% The blue eye is better than the red one... it is the best eye!
