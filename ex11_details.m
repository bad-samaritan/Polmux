% Optical simulation
%
% Example: Details
%

clear all
clc

%addpath('directory where you placed this simulator');   % PATH of .m files

%%%%%%%%%%%%%%%%% Field parameters

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 32; % number of symbols
Nt    = 32; % points x symbol
Nch   = 1;  % number of channels

%%%%%%%%%%%%%%%% Pulse parameters

exratio    = 13;     % extinction ratio [dB]
lam      = 1550;   % central wavelength [nm]
spac     = 0.4;    % channel spacing [nm]
symbrate = 10;     % baudrate [Gbaud]
duty     = 1;      % duty cycle
roll     = 0.2;    % pulse roll-off


%%%%%%%%%% Tx side
outdir = 'pippo';

reset_all(Nsymb,Nt,Nch,outdir);

Pavg = 1;
E    = lasersource(Pavg, lam, spac);

for ii=1:Nch
    pat(:,ii)=pattern('debruijn',ii); % note the different prbs seeds
    elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
    Eopt(:,ii)=mz_modulator(E(:,ii), elec(:,ii),struct('exratio',exratio));
end

create_field('unique',Eopt);


subname = 'Tx'; % suffix for the name
printfield('x',1,subname,'pa--') % print the Tx field (normalized)
                                 % use 'pa--' if you don't want
                                 % normalization

fprintfmsg('Personal message in simul_out');    % Look to simul_out

% Now compress the printed file. Please, check that you have you installed 
% gzip.

system(['gzip ',outdir,'/',outdir,...
    '.MOD/tempx_ch01_',subname,'.dat']); 

% Now read the compressed file. There is no need to uncompress it!

plotfield('x',1,'papa','r-') % plot the Tx field 

% Uncompress the printed file
system(['gunzip ',outdir,'/',outdir,...
    '.MOD/tempx_ch01_Tx.dat']); 

% look at the documentation in html form
mdoc plotfield
