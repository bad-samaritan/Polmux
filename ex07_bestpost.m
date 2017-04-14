% Optical OOK simulation
%
% Example: propagation of a single channel OOK signal within a 20x100 km 
%          of a fully compensated system. Search the best post 
%          compensating fiber that maximizes the eye opening.
%


clear all
clc
%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 32; % number of symbols
Nt    = 32; % points x symbol
Nch   = 1;  % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi      = 0.5*pi; % average cumulated nonlinear phase
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

%%%%%%%%%%%%%%%%% Receiver

x.oftype = 'gauss';   % optical filter type
x.obw    = 1.8;       % optical filter bandwidth 
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = 0.65;      % electrical filter bandwidth
                      % all bandwidths are normalized to the bitrate

%%% Post-compensating fiber parameters
x.slopez = 0;     % post-fiber cumulated slope [ps/nm^2]
x.dpost = [0 1500];  % range for the search
x.lambda = lam;   % post-fiber central wavelength [nm]
x.rec    = 'ook'; % receiver type

% Karhunen-Loeve parameters

x.ber    = 1e-5;       % target bit error rate at which the optimal
                       % post dispersion is searched.
x.eta    = 1.4;        % frequency expansion fatcor (see BER_KL)
x.mu     = 3.5;        % time expansion fatcor (see BER_KL)
x.osnr   = 12+(-1:10); % signal-to-noise ratios [dB/0.1nm]. 
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
Nmid        = ceil(Nch/2);

%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

Pavg = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan); % get the average power
E    = lasersource(Pavg, lam, spac);

for ii=1:Nch % this cycle is not necessary, but it is a good idea to leave
             % it as is, for simple upgrade to WDM
    pat(:,ii)=pattern('debruijn',ii); 
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii)=mz_modulator(E(:,ii), elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt,[],struct('power','average'));

figure(1)
plotfield('x',1,'pa--','r-') % plot the Tx field

%%%%%%%%%% Optical link
for k=1:Nspan
	fiber(tx,'g-sx')
	fiber(comp,'g-sx')
	
	ampliflat(Gerbio,'gain')
end
plotfield('x',1,'pa--','b-') % plot the Rx field without the post-compensation

%%%%%%%%%% Receiver

%%%%%% Now search for the optimal post compensation

[be,bpost] = best_sp(Nmid,x,pat(:,Nmid));

fprintf('\nBest OSNR penalty = %f [dB]\n',be);
fprintf('Best post cumulated dispersion = %f [ps/nm]\n\n',bpost);

% We have a negative penalty of -0.22 dB. Let us look if the result makes
% sense.

Los = 50; % we test 50 post compensating fibers
allpost=linspace(x.dpost(1),x.dpost(2),Los);
x.dpost = 0; % instead of remove if, we fix to 0
post = struct('length',1e3,'alphadB',0,'aeff',80,'n2',2.7e-20,'lambda',...
    lam,'dzmax',1e3,'dphimax',0.1,'slope',0); % we will use an ideal post, 
% hence aeff,n2,dphimax are useless

osnr_loc = zeros(1,Los);
for k=1:Los
    post.disp = allpost(k);
    fiber(post,'g---'); % APPLY POST
    [tp,locos] = ber_kl(Nmid,x,pat(:,Nmid));
    osnr_loc(k) = locos;
    post.disp = -allpost(k);
    fiber(post,'g---'); % REMOVE POST    
end
   
figure(2)
plot(allpost,osnr_loc,'b-o')
grid on; hold on;
xlabel('post cum. dispersion [ps/nm]')
ylabel('OSNR  [dB/0.1nm]')

% Good job? Note how many lines are needed instead of best_sp...
