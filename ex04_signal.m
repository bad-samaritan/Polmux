% Optical simulation
%
% Example: Understanding the signal generation
%

clear all
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64;     % number of symbols
Nt    = 128;    % points x symbol
Nch   = 5;      % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

phi      = 0.2*pi; % peak cumulated nonlinear phase
exratio    = 13;   % extinction ratio [dB]
lam      = 1550;   % central wavelength [nm]
spac     = 0.8;    % channel spacing [nm]
symbrate = 10;     % baudrate [Gbaud]
duty     = 1;      % duty cycle
roll     = 0.2;    % pulse roll-off

% Note: the spectral efficiency is 0.2

%%%%%%%%%%%%%%%% Receiver parameters

oftype = 'gauss';  % optical filter type
obw    = 2;        % optical filter bandwidth 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Tx side

reset_all(Nsymb,Nt,Nch);

Ppeak = 10*(1:Nch);
E    = lasersource(Ppeak, lam, spac);

for ii=1:Nch
    pat(:,ii)=pattern('debruijn',ii);    % note the different de Bruijn seeds
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii)=mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt);

for k=1:Nch
    figure
    plotfield('x',k,'p---','r-',oftype,obw);   % Note the temporary 
                                               % extraction of  channel k.
end

figure  % plot the unique field in frequency. Note the different heights
plotfield('x',1,'--p-')

% you can see different Debruijn sequences. Note that the average power is 1.
% Now, try different manual patterns:

reset_all(Nsymb,Nt,Nch);

Ppeak = 10*(1:Nch);
E    = lasersource(Ppeak, lam, spac);

pat(:,1)=pattern('1100'); 
pat(:,2)=pattern('100'); 
pat(:,3)=pattern('11111111');
pat(:,4)=pattern('1010'); 
pat(:,5)=pattern('1000001'); 

for ii=1:Nch
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii)=mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt);

% Now plot channel by channel
for k=1:Nch
    figure
    hold on;
    plotfield('x',k,'p---','b-',oftype,obw);    % Note the temporary  
                                                % extraction of channel k.
end
% Note that a constant signal is not exactly constant in plotfield, since
% it is extracted with a filter.
