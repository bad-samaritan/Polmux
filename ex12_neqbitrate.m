% Optical simulation
%
% Example: Create signals with different bitrates.
%
% In this example we combine a 40 Gb/s and a 10 Gb/s signal using a trick.
% In principle, such combination is not possible since the global variable
% GSTATE.SYMBRATE is a scalar, but we can exploit the potentiality of the
% code to get our target. 
%

clear all
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 32;  % number of symbols
Nt    = 256; % points x symbol
Nch   = 2;   % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

exratio     = 13;   % extinction ratio [dB]
lam       = 1550;   % central wavelength [nm]
spac      = 2;      % channel spacing [nm]
symbrate1 = 10;     % baudrate 1 [Gbaud]
symbrate2 = 40;     % baudrate 2 [Gbaud]
duty      = 1;      % duty cycle
roll      = 0.2;    % pulse roll-off

%%%%%%%%%%%%%%%% Receiver parameters

oftype = 'butt6'; % optical filter type
obw    = [2 8];   % optical filter bandwidths referred to bitrate OF CHANNEL 1
                  % Remember that after the final generation of the WDM
                  % signal, GSTATE.SYMBOLRATE is set to the last value, here
                  % the bitrate 10 Gb/s.
                  % It turns out that we are using 20 GHz for channel 1
                  % and 80 GHz for channel 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Tx side

% We first generate the 40 Gb/s signal. A 40 Gb/s signal corresponds to a
% 10 Gb/s signal with time compressed of a factor 4. Hence the 40 Gb/s can
% carry 4 times the number of bit of the 10 Gb/s over the same temporal 
% interval. But all signals must have the same length in Optilux. Hence, 
% for instance, we can allocate 4 times the number of bits by dividing Nt  
% by 4. As a result the 40G signals has 4*Nbit bits with Nt/4 points x bit.

reset_all(Nsymb*4,Nt/4,Nch);

for ii=1:Nch
    pat(:,ii)=pattern('debruijn',ii);    
    elec2(:,ii)=electricsource(pat(:,ii),'ook',symbrate2,'cosroll',duty,roll); 
end

% Remember that ELECTRICSOURCE do not use directly the variable bitrate.

% Now we reset_all again and we create two signals at 10 Gb/s.
reset_all(Nsymb,Nt,Nch);

for ii=1:Nch
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate1,'cosroll',duty,roll); 
end

% We now discard the second signal and use the 40 Gb/s signal!!
elec(:,2) = elec2(:,2);

Pavg=1;
E    = lasersource(Pavg, lam, spac);

for ii=1:Nch
  Eopt(:,ii) = mz_modulator(E(:,ii), elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt);

% Let us look to the spectrum
figure(1)
plotfield('x',1,'--p-','b-');

% It is ok.

% Let us now extract the two channels
for k=1:Nch
    figure(k+1)
    plotfield('x',k,'p-p-','b-',oftype,obw(k)); % Note the temporary 
                                                % extraction of  channel k.
end
