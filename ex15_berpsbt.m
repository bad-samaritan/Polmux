% Optical PSBT simulation
%
% Example: Evaluate the bit error rate (BER) of a PSBT signal
%

clear all;
clc;

%%%%%%%%%%%%%%%%% Field parameters

Nsymb = 64; % number of symbols
Nt    = 32; % points x symbol
Nch   = 1;  % number of channels

%%%%%%%%%%%%%%%%  Pulse parameters

exratio    = 13;   % extinction ratio [dB]
lam      = 1550; % central wavelength [nm]
spac     = 0.4;  % channel spacing [nm]
symbrate = 10;   % baudrate [Gbaud]
duty     = 1;    % duty cycle
roll     = 0.2;  % pulse roll-off

%%%%%%%%%%%%%%%%% Receiver

x.oftype = 'gauss'; % optical filter type
x.obw    = 1.5;         % optical filter bandwidth
x.eftype = 'bessel5';    % electrical filter type
x.ebw    = inf;          % electrical filter bandwidth
                         % all bandwidths are normalized to the bitrate
x.rec    = 'ook';        % receiver type
x.plot   = 'ploteye';    % type of plot
x.color  = 'b-';         % color of plot

% Karhunen-Loeve parameters

x.ber    = 1e-9;      % target bit error rate at which the
                      % OSNR is evaluated.
x.eta    = 3;         % frequency expansion fatcor (see BER_KL)
x.mu     = 3.5;       % time expansion factor (see BER_KL)
x.osnr   = 3+(-1:25); % signal-to-noise ratios [dB]
x.poln   = 2;         % noise polarizations        
x.saddle = 'y';       % evaluate the BER without saddlepoint 
                      % approximation (see BER_KL)

reset_all(Nsymb,Nt,Nch);

Pavg = 1;
E    = lasersource( Pavg, lam, spac);

for ii = 1:Nch

    pat(:,ii)    = pattern('debruijn',ii);
    pat_rx(:,ii) = pat_decoder(pat(:,ii),'psbt');

    elec(:,ii)  = electricsource(pat(:,ii),'psbt',symbrate,'cosroll',duty,roll);

    Eopt(:,ii) = mz_modulator(E(:,ii), elec(:,ii));
    
    % Filter PSBT signal to improve OSNR
    % Eopt(:,ii) = lpfilter(Eopt(:,ii),'gauss',0.7);
    % Alternatively you can implement EPSBT
    % elec2(:,ii) = electricsource(pat_rx(:,ii),'ook',symbrate,'cosroll',duty,roll);
    % Eopt(:,ii)  = mz_modulator(Eopt(:,ii),elec2(:,ii),struct('exratio',0));
end

create_field('unique',Eopt);

figure(1)
hold on
plotfield('x',1,'pap-','b',x.oftype,x.obw);

figure(2)
hold on
[pb,osnr]=ber_kl(1,x,pat_rx(:,1));


figure(3)
hold on
plot(x.osnr,log10(pb),'b');
grid on;
xlabel('OSNR   [dB/0.1nm]')
ylabel('BER')
axis([2 15 -15 0])

