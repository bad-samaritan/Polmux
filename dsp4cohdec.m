function [Phases Amplitudes worsteyeop] = dsp4cohdec( ich, pat, x, p)

%DSP4COHDEC Digital signal processing for a coherent receiver
%   [PHASES AMPLITUDES]=DSP4COHDEC(ICH,PAT,X,P) returns the received phase
%   (PHASES) and amplitude (AMPLITUDES) of a quadrature phase shift keying
%   (QPSK) signal after digital signal processing (DSP). 
%   ICH is the channel number, PAT the pattern of symbols. 
%
%   The processing consists in many steps:
%
%   - Analog to Digital Conversion ( quantization )
%   - Dispersion Compensation through FIR
%   - Rebuilding of Complex Signal(s)
%   - Digital Clock Recovery
%   - Polarization Demultiplexing [1-2] 
%   - Frequency Estimation using the Bell Labs alghorithm [1]
%   - Phase Estimation using the Viterbi&Viterbi alghorithm [2]
%   - Hard Decision
%
%   X a struct whose fields are:
%
%   X.oftype = optical filter type (see MYFILTER)
%   X.obw = optical filter 3 dB bandwidth normalized to the symbol rate
%   X.oord = optical filter order (if X.oftype is 'supergaussian')
%   X.eftype = electrical filter type (see MYFILTER)
%   X.ebw = electrical filter 3-dB bandwidth normalized to the symbol rate
%   X.eord = electrical filter order (if X.eftype is 'supergaussian')
%
%   X can also have the optional parameters:
%
%   X.lodetuning    = Local Oscillator Detuning frequency [Hz]
%   X.lophasenoise  = Local Oscillator Phase Noise Vector [rad]
%   X.lolinewidth   = Local Oscillator Linewidth (Hz/symbolrate)
%   X.lopower       = Local Oscillator Power [dBm]
%   X.ts = Fixed sampling time (-0.5 <= X.ts <= 0.5) for the eye. 
%   X.plot = 'ploteye': plots the phase-eye in the active figure; 
%       'plotcur' plots the received current.
%   X.pol = 'x','y','xy'. Polarization of X.plot. In absence of this flag
%       all active polarizations will be plotted.
%   X.color = color string for the plot (see PLOT). E.g. 'b-'.
%   X.dpost = post compensating fiber cumulated dispersion [ps/nm]
%   X.slopez = post compensating fiber cumulated slope  [ps/nm^2]
%   X.lambda = wavelength [nm] at which the post compensating fiber
%              has a cumulated dispersion equal to X.dpost.
%   X.delay = 'theory' means that the delay uses the theoretical delay 
%           saved within GSTATE.DELAY (see CREATE_FIELD). By default the 
%           delay is measured by a cross-correlation measurement between 
%           the received current and an artificial digital signal with 
%           ideal non-return to zero bits with symbols equal to PAT. 
%           The correlation method is useful in presence of polarization 
%           mode dispersion (PMD).
%
%   P is a structure whose fields are:
%
%       P.baudrate      : Baud-Rate or symbol rate [GBaud]
%       P.workatbaudrate: work with one (true) or two (false) samples/symb.
%       P.sps           : samples per symbol in the IRX matrix
%       P.applyadc      : enable quantization (true/false) 
%       P.adcbits       : number of bits of the ADC
%       P.applydcf      : enable dispersion compensation filter (true/false)
%       P.dispersion    : dispersion [ps/nm]
%       P.lambda        : wavelength of the channel [nm]
%       P.ndispsym      : dispersion compensating filter length [symbols]
%       P.applynlr      : enable nonlinear rotation (true/false)
%       P.nlralpha      : constant for non linear rotation algorithm
%       P.applypol      : enable polarization demux/tracking (true/false)
%       P.polmethod     : algorithm: 'singlepol', 'CMA','EASI','Combo'.
%       P.cmaparams.R   : radius of CMA 
%       P.cmaparams.mu  : convergence parameter of CMA
%       P.cmaparams.taps: CMA filter taps
%       P.cmaparams.txpolars: transmitted polarizations
%       P.cmaparams.phizero: rotation angle between tx and rx field 
%       P.easiparams.mu: convergence parameter of EASI
%       P.easiparams.txpolars: transmitted polarizations 
%       P.easiparams.phizero: rotation angle between tx and rx field
%       P.modorder      : modulation order: QPSK=2
%       P.freqavg       : frequency estimation smoothing parameter
%       P.phasavg       : phase estimation smoothing parameter
%       P.poworder      : modulus power in phase estimation. See [4]
%
%   Given PHASES and AMPLITUDES the symbol constellation can be plotted
%   using the function POLAR:
%
%   POLAR(PHASES,AMPLITUDES)
%
%   References:
%   [1] Satoshi Tsukamoto, Yuta Ishikawa, and Kazuro Kikuchi , "Optical
%   Homodyne Receiver Comprising Phase and Polarization Diversities with 
%   Digital Signal Processing", Proceedings of ECOC Conference, Paper
%   Mo4.2.1, ECOC 2006, Cannes, France
%   [2] Dominique N. Godard, "Self-Recovering Equalizationand Carrier
%   Tracking in Two-Dimensional Data Communication Systems", IEEE Trans. on
%   COMMUNICATIONS, VOL. COM-28, No. 11, Nov 1980, pp. 1867-1875
%   [3] A. Leven et al. "Frequency Estimation in Intradyne Reception," IEEE 
%	Photonics Technology Letter, vol. 19, no. 6,pp. 366-368, Mar. 2007  
%   [4] A. J. Viterbi and A. M. Viterbi, "Nonlinear Estimation of 
%	PSK-Modulated Carrier Phase with Application to Burst Digital 
%	Transmission," IEEE Transactions  on Information Theory, vol. IT-29, 
%	no. 4, pp. 543-551, Jul. 1983
%
%   See also PATTERN, MYFILTER, RECEIVER_COHMIX, BER_ESTIMATE
%
%   Author: Massimiliano Salsi, 2009
%   Modified by Marco Bertolini, 2009
%   Modified by Paolo Serena, 2009
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2009  Paolo Serena, <serena@tlc.unipr.it>
%			 
%    Optilux is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    Optilux is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

global GSTATE
global CONSTANTS;  % CONSTANTS is a global structure variable.
CLIGHT = CONSTANTS.CLIGHT;  % speed of light [m/s]

% Photodiodes Order:
% Channel 1 : Polar 1, In phase
% Channel 2 : Polar 2, In phase
% Channel 3 : Polar 1, In quadrature
% Channel 4 : Polar 2, In quadrature
COS_POL1 = 1;
COS_POL2 = 3;
SIN_POL1 = 2;
SIN_POL2 = 4;

switch x.rec
    case 'coherent'        
        [Irxt,x] = receiver_cohmix(ich,x); % get the current
        % if GSTATE.FIELDY is created along propagation, it doesn't carry
        % information, thus we must ignore it and consider only the X
        % polarization for which the pattern is defined
    otherwise
        error('Flag X.rec must be ''coherent''');
end
if isempty(GSTATE.FIELDY) || (size(pat,2) == 1)
    Irx = Irxt(:,1:2); % use X polarization
    x.isy = 0;
else
    Irx = Irxt;        % use both polarizations
    x.isy = 1;
end

% ADC finite bit number emulator
if p.applyadc
    M = max( max( abs( Irx ) ) );   % <- Optimal range. There's no loss.
    % Quantization:
    Irx = round( ( Irx + M ) / 2 ./ M * 2^p.adcbits )* 2 .* M / 2^p.adcbits - M;
end
[eyeb,best_ts,delay]=mygeteyeinfo(ich,Irx,x,pat);
for npol = 1:size(Irx,2)/2
    Irx(:,(1:2)+2*(npol-1)) = fastshift(Irx(:,(1:2)+2*(npol-1)),round(-delay(npol)*GSTATE.NT));
end
% we will sample in the middle of the period

% The ADC works at a given sampling rate. This aspect cannot be neglected
% like the finite bit number. So now there is the mandatory decimation:
if p.workatbaudrate
    DecimationRate = p.sps; % We will have one sample x symbol
else
    DecimationRate = p.sps / 2; % We will have two samples x symbol
end
if size( Irx,2) == 2
    Irxdec = [ decimate( Irx(:,1), DecimationRate,16,'fir' ) ...
            decimate( Irx(:,2), DecimationRate,16,'fir' ) ];
else
    Irxdec = [ decimate( Irx(:,1), DecimationRate,16,'fir' ) ...
            decimate( Irx(:,2), DecimationRate,16,'fir' ) ...
            decimate( Irx(:,3), DecimationRate,16,'fir' ) ...
            decimate( Irx(:,4), DecimationRate,16,'fir' ) ];
end
% Rebuilding Signals from Sequences ____________________________________
switch size( Irxdec, 2 ) 
    case 2
        Signals = complex( Irxdec(:,1), Irxdec(:,2) );
    case 4
        Signals = [ complex( Irxdec(:,COS_POL1), Irxdec(:,SIN_POL1) ) ...
                    complex( Irxdec(:,COS_POL2), Irxdec(:,SIN_POL2) ) ];
    otherwise
        error('Unable to manage Irx matrix.');
end
% _____________________________________________________________________

% Digital Dispersion Compensation______________________________________
if p.applydcf
    Beta2L           = -p.dispersion * p.lambda^2 / 2 / pi / CLIGHT * 1E-21;
    SamplesPerSymbol = 1 + ~p.workatbaudrate;
    DCF_Points       = length(Signals);
    DCF_Length       = p.ndispsym*SamplesPerSymbol;
    DCF_BandW        = SamplesPerSymbol*p.baudrate; % bitrate in Hz
    Hfilt            = DispCompFilter( Beta2L, DCF_BandW, DCF_Points , DCF_Length );
    Hfilt = Hfilt * ones(1,size(Signals,2));
    Signals     = ifft( fft(Signals) .* Hfilt );
end
% ___________________________________end of Dispersion Compensation

% DIGITAL CLOCK RECOVERY_______________________________________________
% This operation is useless for the moment. We are already sampling in the
% middle of the period. Now that we have compensated for chromatic
% dispersion we need only the central samples.
if ~p.workatbaudrate
    Signals = Signals(1:2:end,:);
end
%______________________________________________________________________

% NON LINEAR PHASE ROTATION____________________________________________
if p.applynlr
    Signals    = NLRotation(Signals, p.nlralpha);
end
%______________________________________________________________________
% NORMALIZATION
peak = 4*sqrt(GSTATE.POWER(ich)); % 2* -> see receiver_cohmix
Signals = Signals/peak;

%______________________________________________________________________
% POLAR ROTATION: None, Fixed or Adaptative ___________________________
if p.applypol && (size(Signals,2)==2)
    switch lower(p.polmethod)
        case 'singlepol'
            r = Signals(:,1) ./ Signals(:, 2);
            rKikuchi = mean( r );
            Signals = rotpolar(Signals, rKikuchi );
        case 'cma'
            Signals = cmapolardemux( Signals, p.cmaparams );
        case 'easi'
            Signals = easipolardemux( Signals, p.easiparams );
        case 'combo'
            Signals = easipolardemux( Signals, p.easiparams );
            Signals = cmapolardemux( Signals, p.cmaparams );
        otherwise
            error('Unknown Polar Rotation method.');
    end
end
%_________________________________________________END OF POLAR ROTATION

%   CARRIER PHASE ESTIMATION_______________________________________________
M = 2.^p.modorder;
navg = p.freqavg;
if navg
    % Estimating frequency in bell labs style [citation neede]:
    omega = cumsum( vitvit( Signals .* conj( fastshift(Signals,1) )...
        , M, M, navg , false) );
    % Cleaning omega to match the circularity:
    closestallowedendpoints = omega(1,:) + round((omega(end,:)-omega(1,:))/2/pi)*2*pi;
    correctionratio = closestallowedendpoints ./ omega(end,:);
    omega=((omega-ones(length(omega),1)*omega(1,:)).*(ones(length(omega),1)*correctionratio))...
        +ones(length(omega),1)*omega(1,:);
    % Demodulating signals:
    sigdemod = Signals .* fastexp( -omega );
    % Estimating phase using Viterbi and Viterbi method:
    navg    = p.phasavg;
    P       = p.poworder;
    theta   = vitvit( sigdemod, P, M, navg, true );
    if p.modorder > 1
        CarrierPhaseOffSet = +pi/4;
    else
        CarrierPhaseOffSet = 0;
    end
    Carrier = fastexp( -omega - theta + CarrierPhaseOffSet );
else
    navg    = p.phasavg;
    P       = p.poworder;
    theta   = vitvit( Signals, P, M, navg, true );
    if p.modorder > 1
        CarrierPhaseOffSet = +pi/4;
    else
        CarrierPhaseOffSet = 0;
    end
    Carrier = fastexp(- theta + CarrierPhaseOffSet );
end
Signals = Signals .*  Carrier ;

Amplitudes  = abs( Signals );
Phases      = angle( Signals );
worsteyeop = min(eyeb(abs(eyeb-mod(eyeb,pi/2))<1e-10));

% pat_rx = extractPattern( Phases );
% 
% switch p.modorder
%     case 1
%         pat_rx = diff_decod(pat_rx);
%     case 2
%         pat_rx = [ diff_decod2(pat_rx(:,1:2)) diff_decod2(pat_rx(:,3:4)) ];
%     otherwise
%         error('dsp4cohdec: error during decoding.');
% end

return ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hfilt = DispCompFilter( Beta2L, B , N, FilterLength)
freq = ( -B/2 : B / N : B/2 * (N-2) / N );
freq = ifftshift( freq );
delay = 2*pi*freq/B * ( FilterLength/2  );
argum = (2*pi*freq).^2 * Beta2L / 2 - delay;
H = fastexp ( argum );
b = ifft( H );
b = b(1:FilterLength + 1 );
Hfilt = fft(b.', N) .* fastexp(delay.');
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = NLRotation( s, alpha)
Phases     = angle( s );
Amplitudes = abs( s );
Asquare    = sum( Amplitudes .* Amplitudes , 2 );
Delta_P    = Asquare - mean( Asquare );
Phases     = Phases + alpha*Delta_P*ones(1,size(s,2));
y          = Amplitudes .* fastexp( Phases );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = vitvit( s, P, M, k, applyunwrap)
L = length(s);
% Phase -> Phase * M
% Amplitude -> Amplitude ^ P
if P == M
    s = s .^ P;
else
    s = abs(s).^P .* fastexp( angle( s.^ M ) );
end
if k>0
N = 2*k + 1;
if N<L
    Smoothing_Filter = fft( ones(N, 1) / N , L ) * ones(1,size(s,2));
    s = ifft( fft(s,L) .* Smoothing_Filter );
else
    slong = repmat( s, ceil(N / L), 1 );
    Smoothing_Filter = fft( ones(N, 1) / N , ceil(N / L).*L ) * ones(1,size(s,2));
    slong = ifft( fft(slong) .* Smoothing_Filter );
    s = slong(1:L,:);
end
end
if applyunwrap
    theta = unwrap( angle(s) ) / M;
else
    theta = angle( s ) / M;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = rotpolar( x , r )
if abs(r) < 0.5
    m = abs(r);
    delta = angle(r);
    alpha = m^2 / (m^2 + 1);
else
    m = abs( 1 / r);
    delta = -angle( 1 / r );
    alpha = 1 / ( m^2 + 1 );
end
M = [sqrt(alpha)*exp(-j*delta)   -sqrt(1-alpha)*exp(-j*delta);...
     sqrt(1-alpha)             sqrt(alpha)   ];
y = x * M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pattern=extractPattern(phase)
% second_bit = phase>0;
% first_bit = abs( phase) <= pi/2;
% if size(phase,2) == 1
%     pattern = [first_bit second_bit];
% else
%     pattern = [first_bit(:,1) second_bit(:,1)  first_bit(:,2)  second_bit(:,2) ];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=cmapolardemux(x, params)
R       = params.R;                 % Radius for convergence
mu      = params.mu;                % Convergence parameter
taps    = params.taps;              % length of adaptive filters
halftaps= floor( taps / 2);          % half of (taps-1)
polmux  = (params.txpolars == 2); % 1 o 2 polarization transmitted
hzero = zeros( taps , 2, 2);
if isfield(params,'mat')
    M = params.mat;
else
    if ~polmux
        r = mean( x(:,1) ./ x(:,2) );
        M = rotpolar(1, r).';
    else
        phi = params.phizero;
        M = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    end
end
hzero( halftaps+1, :, :) = M;  % Initializing central taps
if halftaps
    extendedx = [ x(end-halftaps+1:end,:) ; x ; x(1:halftaps,:) ];
else
    extendedx = x;
end
h1 = squeeze( hzero(:, 1, :) );
h2 = squeeze( hzero(:, 2, :) );
if taps==1
    h1 = h1.';
    h2 = h2.';
end
convergence = false;
c=1;
L=length(x);
repetitions = 50*ceil(1./(L.*mu));
while ~convergence && (c<(repetitions))

    h1_old = h1*1;
    h2_old = h2*1;
    [y h1_new h2_new] = cmaadaptivefilter(extendedx, h1, h2, taps, mu, R, 1);
    % Here I'm fixing to 1 the number of samples per bit. This can be 2 or
    % more, but syill is not immplemented in dsp4cohdec. In case this
    % changes, the 1 must be replaced with the correct value.
    if any(any(h1_new)) || any(any(h2_new))
        h1=h1_new;
        h2=h2_new;
    end
    if max(max(abs([h1_old-h1 h2_old-h2]))) < 5e-5 % max(max(abs([(h1_old-h1)./(h1+h1_old)/2 (h2_old-h2)./(h1+h1_old)/2])))
        convergence = true;
    end
    c=c+1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=easipolardemux(x, params)
mu      = params.mu;                % Convergence parameter
taps    = 1;             
halftaps= floor( taps / 2);          % half of (taps-1)
polmux  = (params.txpolars == 2); % 1 o 2 polarization transmitted
hzero = zeros( taps , 2, 2);
if isfield(params,'mat')
    M = params.mat;
else
    if ~polmux
        r = mean( x(:,1) ./ x(:,2) );
        M = rotpolar(1, r).';
    else
        phi = params.phizero;
        M = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    end
end
hzero( halftaps+1, :, :) = M;  % Initializing central taps
if halftaps
    extendedx = [ x(end-halftaps+1:end,:) ; x ; x(1:halftaps,:) ];
else
    extendedx = x;
end
h1 = squeeze( hzero(:, 1, :) );
h2 = squeeze( hzero(:, 2, :) );
if taps==1
    h1 = h1.';
    h2 = h2.';
end
convergence = false;
c=1;
L=length(x);
repetitions = 20*ceil(1./(L.*mu));
while ~convergence && (c<(repetitions))

    h1_old = h1*1;
    h2_old = h2*1;
    [y h1_new h2_new] = easiadaptivefilter(extendedx, h1, h2, taps, mu, 1);
    % Here I'm fixing to 1 the number of samples per bit. This can be 2 or
    % more, but syill is not immplemented in dsp4cohdec. In case this
    % changes, the 1 must be replaced with the correct value.
    if any(any(h1_new)) || any(any(h2_new))
        h1=h1_new;
        h2=h2_new;
    end
    if max(max(abs([h1_old-h1 h2_old-h2]))) < 5e-5 % max(max(abs([(h1_old-h1)./(h1+h1_old)/2 (h2_old-h2)./(h1+h1_old)/2])))
        convergence = true;
    end
    c=c+1;

end
function [eyeb,best_ts,delay,xopt]=mygeteyeinfo(ich,Iric,x,pat)

% Eye evaluations. On input:
%   ICH: channel number
%   IRIC: electric current returned by the receiver
%   X: struct variable (see EVAL_EYE)
%   PAT: symbols pattern
%
%   On output:
%
%   EYEB: Eye opening
%   BEST_TS: best normalized sampling time (center of bit is 0)
%   DELAY: overal normalized (to the bit time) delay
%   XOPT: Best sampling time in discrete points.
%

global GSTATE


for npol = 1:1+x.isy
    ipat = pat(:,npol);

    ipat_old = ipat;
    ipat(ipat == 0) = -3/4*pi;
    ipat(ipat == 1) = 3/4*pi;
    ipat(ipat == 2) = -1/4*pi;
    ipat(ipat == 3) = 1/4*pi;

    %         Nfft = GSTATE.NSYMB*GSTATE.NT;
    %         refsig = reshape(repmat(ipat,1,GSTATE.NT)',Nfft,1);
    if isfield(x,'delay') && strcmp(x.delay,'theory') % theoretical delay
        if isfield(x,'b2b')
            if strcmp(x.b2b,'b2b'), avgdelay = 0;end   % back-to-back: no delay
        else
            if x.isy
                avgdelay = 0.5*(GSTATE.DELAY(1,ich)+GSTATE.DELAY(2,ich));
            else
                avgdelay = GSTATE.DELAY(1,ich);
            end
        end
        delay = ones(1,1+x.isy)*(avgdelay+evaldelay(x.oftype,x.obw*0.5)+...
            evaldelay(x.eftype,x.ebw)+x.post_delay);
        Iric_t = angle(complex(Iric(:,1+2*(npol-1)),...
            Iric(:,2+2*(npol-1))));
    else
        [delay(npol),wrn,rho,Iric_t] = corrdelay(complex(Iric(:,1+2*(npol-1)),...
            Iric(:,2+2*(npol-1))),ipat,GSTATE.NT,GSTATE.NSYMB,'phase');
    end
    halfbit = GSTATE.NT/2;

    nshift = round(halfbit-delay(npol)*GSTATE.NT); % the first bit is centered at index 1
    Iricmat = reshape(fastshift(Iric_t,nshift),GSTATE.NT,GSTATE.NSYMB*size(Iric_t,2))';  % Note the transpose!

    for nalph = 1:max(ipat_old)+1
        minv(:,nalph+(max(ipat_old)+1)*(npol-1)) = min(Iricmat(ipat_old==nalph-1,:),[],1);
        maxv(:,nalph+(max(ipat_old)+1)*(npol-1)) = max(Iricmat(ipat_old==nalph-1,:),[],1);
    end

    if npol == 1
        Iricmatx = Iricmat;
        Iricx       = Iric_t;
    else
        Iricmaty = Iricmat;
        Iricy       = Iric_t;
    end
end

if (~isempty(minv) && ~isempty(maxv))

    eyeop(:,1) = minv(:,3)-maxv(:,1);    % eye opening
    eyeop(:,2) = minv(:,2)-maxv(:,4);
    eyeop(:,3) = minv(:,4)-maxv(:,3);
    eyeop(:,4) = minv(:,1)-(maxv(:,2)-2*pi);
    if x.isy
        eyeop(:,5) = minv(:,7)-maxv(:,5);    % eye opening
        eyeop(:,6) = minv(:,6)-maxv(:,8);
        eyeop(:,7) = minv(:,8)-maxv(:,7);
        eyeop(:,8) = minv(:,5)-(maxv(:,6)-2*pi);
    end
    worsteye = min(eyeop,[],2);

    if isfield(x,'ts') % fixed sampling time
        xopt = round((x.ts+0.5)*GSTATE.NT);
        eyeb = eyeop(xopt,:);
        eyeb(eyeb<0) = NaN;
        best_ts = x.ts;
    else
        [eyeb,best_tsn] = max(worsteye);
        eyeb =eyeop(best_tsn,:);

        if GSTATE.NT == 2
            best_ts = best_tsn/GSTATE.NT-0.5;
        else    % Interpolation for the best eye opening
            if best_tsn == 1
                to = 1:3;
            elseif best_tsn == GSTATE.NT
                to = best_tsn-2:best_tsn;
            else
                to = best_tsn-1:best_tsn+1;
            end
            % Parabolic interpolation around the max using three points.
            xopt = 0.5*(to(3)^2*(worsteye(to(1))-worsteye(to(2)))+to(1)^2*(worsteye(to(2))-worsteye(to(3)))+...
                to(2)^2*(worsteye(to(3))-worsteye(to(1))))/(to(3)*(worsteye(to(1))-worsteye(to(2)))+...
                to(1)*(worsteye(to(2))-worsteye(to(3)))+to(2)*(worsteye(to(3))-worsteye(to(1))));
            best_ts = xopt/GSTATE.NT-0.5;
            for nalph = 1:size(eyeop,2)
                eyeb(nalph) = (xopt-to(2))*(xopt-to(3))/((to(1)-to(2))*(to(1)-to(3)))*eyeop(to(1),nalph)+...
                    (xopt-to(1))*(xopt-to(3))/((to(2)-to(1))*(to(2)-to(3)))*eyeop(to(2),nalph)+...
                    (xopt-to(1))*(xopt-to(2))/((to(3)-to(1))*(to(3)-to(2)))*eyeop(to(3),nalph);
            end
        end
    end
else     % if (~isempty(min1) && ~isempty(max0))
    eyeb = Inf;     % the eye does not exist
    best_ts = 0;    % any time is good
    xopt = 0.5*GSTATE.NT;
end


%%%%% GRAPHICAL OUTPUT

if isfield(x,'plot')
    if ~isfield(x,'color')
        col = 'b-';
    else
        col = x.color;
    end

    if  strcmpi(x.plot,'ploteye')
        ntime = -0.5:1/(GSTATE.NT-1):0.5;
        if ~isfield(x,'pol')
            if x.isy, 
                x.pol = 'xy';
            else
                x.pol = 'x';
            end
        end
        if isfield(x,'pol')
            if strcmp(x.pol,'y')
                plot(ntime,Iricmaty',col)
            elseif strcmp(x.pol,'xy')
                subplot(1,2,1)
                plot(ntime,Iricmatx',col);
                subplot(1,2,2)
                plot(ntime,Iricmaty',col)
            else
                plot(ntime,Iricmatx',col)
            end
        end
        xlabel('time')
        ylabel('Eye')
    elseif strcmpi(x.plot,'plotcur')
        time=0:1/GSTATE.NT:GSTATE.NSYMB-1/GSTATE.NT;
        if isfield(x,'pol') && strcmp(x.pol,'y')
            plot(time,Iricy,col);
        else
            plot(time,Iricx,col);
        end
        xlabel('time   [a.u.]');
        ylabel('Current   [mA]')
    elseif strcmp(x.plot,'')
        % nothing to do
    else
        error(['the field plot must be ',...
            '''ploteye'' or ''plotcur''']);
    end
    drawnow;
end
