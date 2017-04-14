% TEST_DOP
% test the Optilux function DOPmeter() with 1 channel
% This example shows how the signal to noise ratio can be evaluated by a
% DOP measure, in the case of totally polarized field affected by
% unpolarized ASE noise.
%
% Armando Vannucci, 2009

clear all
clc

%%%%%%%%%%%%%%%%% Field parameters
Nsymb = 32;     % number of bits
Nt = 64;        % points x bit
Nch = 1;        % number of channels
Npol = 2;       % number of polarizations
%%%%%%%%%%%%%%%%  parameters continue below

reset_all(Nsymb,Nt,Nch); 
global CONSTANTS;            % CONSTANTS is a global structure variable.

%%%%%%%%%%%%%%%%  Pulse parameters
symbrate = 10;   % bit rate [Gb/s]
duty = 1;        % duty cycle
roll = 0.2;      % pulse roll-off
exratio = 12;    % extinction ratio [dB]
lam = 1550;      % central wavelength [nm]
spac=0.4;        % channel spacing [nm]
P_S= 1.0;        % [mW] average signal power
%%%%%%%%%%%%%%%%  Opt.Ampli parameters
ampgain= 20;     % [dB] dummy gain: used only for generating noise
Glin= 10^(ampgain/10); % linear gain
% uncomment this for a large DOP (and a small cloud of pol.states)
P_ASE= 0.001;              % ASE noise variance (power [mW])
% uncomment this for a DOP=0.5 (and a cloud of pol.states)
% P_ASE= 1.0;              % ASE noise variance (power [mW])
% P_S will be amplified by Glin, so we amplify P_ASE to get OSNR= P_S/P_ASE on the total simulation bandwidth
ampP_ASE= P_ASE*Glin;
% Noise factor: F= P_ASE/(hu Btot*(G-1)) 
huBtot= CONSTANTS.HPLANCK*(CONSTANTS.CLIGHT/lam)*(Nt*symbrate)*1e21;   
% *1e21, since lam [nm] and symbrate [Gs/s] and huBtot [mW]
ampopt.f = 10*log10( ampP_ASE/(huBtot*(Glin-1)) ); % noise figure [dB]
% the one-sided ASE PSD of a single amplifier is 2*nsp*hvdl*(Gerbio-1) 
randn('state',3);   % initialize random gaussian generator
%%%%%%%%%%%%%%%%  end parameters

pat=pattern('debruijn',1);
elec= electricsource(pat,'ook',symbrate,'cosroll',duty,roll);


E = lasersource(P_S, lam, spac);
Eoptx = mz_modulator(E,elec,struct('exratio',exratio));
Eopty = zeros(size(Eoptx));

chi= pi/8; phi= 3*pi/4;     % pol. state of polarized CW (linear diagonal -45deg)
[Eoptx,Eopty] = set_sop(Eoptx,Eopty,chi,phi,'aarphd'); % set the SOP
create_field('unique',Eoptx,Eopty,struct('power','average'));


figure(1)
plotfield('tot',1,'p---','r-')    % plot the Tx field (normalized)
                                  % use 'n---' if you don't want
                                  % normalization
                               
% ADD NOISE: two amplifiers with opposite gains=> no gain; only noise
ampliflat( ampgain,'gain',ampopt);
ampliflat(-ampgain,'gain');

figure(2)
plotfield('tot',1,'p---','r-')    % plot the Tx field (normalized)

myDOP= dop_meter(1);
% NOTE: in theory (if the number of samples is very large), the power of
% the polarized component (P_S) should be a fraction P_S/P_avg of the
% sphere radius P_avg= P_S+P_ASE (unpolarized component)

fprintf('DOP= %f \n', myDOP);
fprintf('Real SNR= %.2f dB; SNR evaluated from DOP= %.2f dB \n',...
    10*log10(P_S./P_ASE), 10*log10(myDOP./(1-myDOP)));
