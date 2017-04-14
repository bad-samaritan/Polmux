% Optical DPSK simulation
%
%  Example: propagation of a single channel OOK signal within a 20x100 km 
%           of a fully compensated fiber. Test the accuracy of best_sp.m 
%           and its problems.
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
exratio    = 13;     % extinction ratio [dB]
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
x.slopez = 0;    % post-fiber cumulated slope [ps/nm^2]
x.lambda = lam;  % post-fiber central wavelength [nm]
x.dpost  = -200; % post-fiber cumulated dispersion [ps/nm]

x.rec = 'ook';                  % receiver type

% Karhunen-Lo�ve parameters

x.ber    = 1e-9;       % target bit error rate at which the optimal
                       % post dispersion is searched.
x.eta    = 1.4;        % frequency expansion fatcor (see ber_kl.m)
x.mu     = 3.5;        % time expansion fatcor (see ber_kl.m)
x.osnr   = 12+(-1:10); % signal-to-noise ratios [dB]. 
x.poln   = 2;          % noise polarizations        
x.saddle = 'y';        % evaluate the BER by saddlepoint 
                       % approximation (see ber_kl.m)

% Note that the OSNR in b2b is 12.1 dB, hence we search the OSNR penalty
% within the range -1:10 around 12 dB. By this way we are able to measure
% penalty up to 9.9 dB

% To gain feeling with x.eta and x.mu try to change them.
% Try also to test the accuracy of the saddlepoint approximation by
% removing the line x.saddle = 'y'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Conversions

gam = 2*pi*tx.n2./(lam.*tx.aeff)*1e21;                            % gamma index [1/W/m]
comp.length = (Din - tx.disp*tx.length/1e3)/comp.disp*1e3;        % comp. fiber length [m]
Gerbio = tx.alphadB*tx.length/1e3 + comp.alphadB*comp.length/1e3; % ampli gain [dB]

%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

Pavg = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan); % get the average power
E    = lasersource(Pavg, lam, spac);

for ii=1:Nch
    pat(:,ii)=pattern('debruijn',ii);
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii) = mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt);

%%%%%%%%%% Optical link
for k=1:Nspan
	fiber(tx,'g-sx')
	fiber(comp,'g-sx')
	
	ampliflat(Gerbio,'gain')
end

%%%%%%%%%% Receiver


%%%%%% now search for the optimal post compensation
x.dpost = [300 1300]; % range for the search

[be,bpost] = best_sp(1,x,pat);

fprintf('\nBest OSNR penalty = %f [dB]\n',be);
fprintf('Best post cumulated dispersion = %f [ps/nm]\n\n',bpost);

% Is it the correct value? Let us answer by exhaustive search of
% the optimal post:

x.dpost = 0; % remove the post inside BER_KL
dpost = 300:10:1300;
be1 = zeros(1,length(dpost));
post = comp;
post.length = 1E3;
post.alphadB = 0;

for k=1:length(dpost)
    post.disp=dpost(k);
    fiber(post,'g---'); % apply post
    [pb1,be1(k)] = ber_kl(1,x,pat);
    post.disp=-dpost(k);
    fiber(post,'g---'); % remove post
end

[be pos]=min(be1);
fprintf(['Best post cumulated dispersion (exhaustive search +/- 10 ps/nm)',...
    '= %f [ps/nm]\n\n'],dpost(pos));
figure(1)
plot(dpost,be1)
grid on;
hold on;

% Ok, the previous answer was right. However the function BEST_SP can fail 
% sometimes. For instance, it fails if the range of the OSNR is too small,
% so that the evaluation of the OSNR penalty returns NaN and the golden
% search algorithm do not know in which direction search the best value.
% The function might also fail when the initial bracket is too large. For
% instance, with the following bracket: 

x.dpost = [-10000 10000];
[be,bpost] = best_sp(1,x,pat);

fprintf('\nBest OSNR penalty = %f [dB]\n',be);
fprintf('Best post cumulated dispersion = %f [ps/nm]\n\n',bpost);

% we have the wrong answer!!! Be careful!
