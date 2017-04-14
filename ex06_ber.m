% Optical OOK simulation
%
% Example: propagation of a single channel OOK signal within 20x100 km of 
%          a fully compensated system. Evaluate the bit error rate (BER).
%


clear all
clc
%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 32; % number of symbols
Nt    = 32; % points x symbol
Nch   = 1;  % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi      = 0.2*pi; % average cumulated nonlinear phase
exratio    = 13;   % extinction ratio [dB]
lam      = 1550;   % central wavelength [nm]
spac     = 0.4;    % channel spacing [nm]
symbrate = 10;     % baudrate [Gbaud]
duty     = 1;      % duty cycle
roll     = 0.2;    % pulse roll-off

%%%%%%%%%%%%%%%% Link parameters

Din   = 0;          % in-line dispersion x span [ps/nm]
Nspan = 20;         % number of spans

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


%%%%%%%%%%%%%%%%% Receiver

x.oftype = 'gauss';   % optical filter type
x.obw    = 1.8;       % optical filter bandwidth 
x.eftype = 'bessel5'; % electrical filter type
x.ebw    = 0.65;      % electrical filter bandwidth
% all bandwidths are normalized to the bitrate

x.rec    = 'ook';     % receiver type
x.plot   = 'ploteye'; % type of plot
x.color  = 'r-';      % color of plot

% Karhunen-Lo�ve parameters

x.ber    = 1e-5;       % target bit error rate at which the
                       % OSNR is evaluated.
x.eta    = 1.4;        % frequency expansion fatcor (see BER_KL)
x.mu     = 3.5;        % time expansion fatcor (see BER_KL)
x.osnr   = 12+(-1:10); % signal-to-noise ratios [dB/0.1nm]
x.poln   = 2;          % noise polarizations        
x.saddle = 'y';        % evaluate the BER by saddlepoint 
                       % approximation (see BER_KL)

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
    pat(:,ii)=pattern('debruijn',ii*2);
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii) = mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt,[],struct('power','average'));

figure(1)
plotfield('x',1,'na--','r-') % plot the Tx field (normalized to Tx peak 
                             % power).

%%%%%%%%%% Optical link
for k=1:Nspan
    fiber(tx,'g-sx')
    fiber(comp,'g-sx')

    ampliflat(Gerbio,'gain')
end
plotfield('x',1,'na--','b-') % plot the Rx field (normalized to Tx peak 
                             % power).

%%%%%%%%%% Receiver

figure(3)   % activate figure(3) for the eye diagram
grid on;
hold on;

% Let us first take a look at the OSNR in back-to-back. Our range for
% x.osnr should be around this value.

x.b2b = 'b2b';
[pb,osnr_b2b]=ber_kl(1,x,pat);

% Now we return to our system in the nonlinear regime.
x=rmfield(x,'b2b');
x.color = 'b-'; % change color for the eye in nonlinear regime

[pb,osnr]=ber_kl(1,x,pat);
title('RED: b2b.   BLUE: after propagation');

fprintf('\n\n=========== Results ===========\n\n');
fprintf('OSNR (b2b) @ BER = %.e : %.4f [dB/0.1 nm]\n',x.ber,osnr_b2b);
fprintf('OSNR with propagation @ BER = %.e : %.4f [dB/0.1 nm]\n',x.ber,osnr);
fprintf('OSNR penalty@ BER = %.e : %.4f [dB]\n\n',x.ber,osnr-osnr_b2b);

figure(4)
semilogy(x.osnr,pb);
grid on;
xlabel('OSNR   [dB/0.1nm]')
ylabel('BER')

% Note that with our choice of x.osnr the BER intersects the horizontal line
% pb=1e-5. Hence the function returned the correct value of the OSNR @
% BER=1e-5, because the interpolation works correctly (see BER_KL). 
% Otherwise, we would had have OSNR @ BER=1e-5 = NaN. Let us demonstrate
% this problem:

x.osnr = 12 + (-1:.2:1);    % Previously, we measured  a penalty of 1.47 dB 
                            % outside this range.

x=rmfield(x,'plot');   % stop to plot the eye
[pb,osnr1]=ber_kl(1,x,pat);

fprintf('OSNR with reduced range @ BER = %.e : %.4f [dB/0.1 nm]\n',...
    x.ber,osnr1);

% And we have OSNR penalty = NaN. The best way to encompass this problem is
% to use a large range for the OSNRs. Alternatively, a faster way uses
% x.extrap='yes', so that the BER can be extrapolated outside the range
% x.osnr, by numerical interpolation (see interp1.m). Let us check this 
% option:

x.extrap = 'yes'; % x.osnr is still outside the range of penalty=1.47 dB
[pb,osnr_ext]=ber_kl(1,x,pat);

fprintf('OSNR with extrap=''yes'' @ BER = %.e : %.4f [dB/0.1 nm]\n',...
    x.ber,osnr_ext);
fprintf('error of extrap: %.4f [dB]\n\n',osnr_ext-osnr);

% However, as in any interpolation scheme, extrapolate a value can be
% dangerous and can give the wrong answer if the function is not well
% approximated by a spline. Be careful!
