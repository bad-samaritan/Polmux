% TEST_polscrambler
% test the Optilux function pol_scrambler('rand') with 1 channel.
% This example shows how the SOP of a single completely polarized signal
% can be randomly moved on the Poincaré sphere by applying polarization
% scrambling.
% The signal is CW, i.e. modulated by a sequence of '1' symbols
% polscrambler('rand',coh_timeR) rotates signal samples separately, in time,
% by randomly changing its parameters every coh_timeR [symbol periods] (coherence time).
% Then, DOP is evaluated, which will eventually approach zero, i.e., complete 
% depolarization, for a large number of samples, since polscrambler()
% uniformly covers the Poincaré sphere.
% NOTE: polscrambler('fixed') can be used to perform a deterministic rotation.
%
% Armando Vannucci, University of Parma, 2009

clear all;


%%%%%%%%%%%%%%%%% Field parameters
Nsymb = 32;      % number of bits
Nt = 32;       % points x bit
Nch = 1;        % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters
symbrate = 10;  % bit rate [Gb/s]
duty = 1;       % duty cycle
roll = 0.2;     % pulse roll-off
exratio = inf;  % extinction ratio [dB]: very small so that zeros and one 
                % lie close to the Poincaré sphere
lam = 1550;     % central wavelength [nm]
spac=0.4;       % channel spacing [nm]
P_S= 1.0;       % [mW] average signal power
%%%%%%%%%%%%%%%%  pol_scrambler parameters
% uncomment this to scramble each time sample independently: output DOP=0
 coh_timeR = 1/Nt;   % scramble every signal sample indipendently 
% uncomment this to scramble each time sample indipendently: output DOP=1
% coh_timeR = Nsymb;   % scramble all signal samples equally
% rand('state',5);   % initialize random gaussian generator for pol_scrambler()
%%%%%%%%%%%%%%%%  end parameters

reset_all(Nsymb,Nt,Nch);    
pat=pattern('debruijn',1);  % all-ones sequence: CW
elec= electricsource(pat,'ook',symbrate,'cosroll',duty,roll);


E = lasersource(P_S, lam, spac);
Eoptx = mz_modulator(E,elec,struct('exratio',exratio));
Eopty = zeros(size(Eoptx));

chi= pi/8; phi= 3*pi/4;     % pol. state of polarized CW (linear diagonal -45deg)
[Eoptx,Eopty] = set_sop(Eoptx,Eopty,chi,phi,'aarphd'); % set the SOP
create_field('unique',Eoptx,Eopty,struct('power','average'));

myDOP= dop_meter(1);
avgE= avg_power(1,'abs');
fprintf('DOP before scrambling= %f . mean power = %f mW \n', myDOP, avgE);

pol_scrambler('rand',coh_timeR);  % scrambles ALL samples on the same random SOP

myDOP= dop_meter(1);
avgE= avg_power(1,'abs');
fprintf('DOP after scrambling= %f . mean power = %f mW \n', myDOP, avgE);

