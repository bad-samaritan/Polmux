function varargout = receiver_dqpsk(ich,x)

%RECEIVER_DPSK Complete DQPSK receiver. (POST fiber+OBPF+MZ+LPF).
%   VARARGOUT=RECEIVER_DPSK(ICH,X) returns the received current of a DQPSK 
%   transmission using the following receiver:                            
%
%
%                                Mach Zehnder     __|__
%                                 ---------        / \
%                                /  -pi/4   \   /  ---     -------
%                               /            \ /    |     |       | IRIC(:,1)                            
%                             -               -     |-----|  LPF  |-------                                                   
%          __                /  \   ------   / \  __|__   |       |   
%         /  \              /    \_| 1 bit|_/   \  / \     -------
%        |    |  -------   /       | delay|        ---
%  sig    \  /  |       | /         ------          |                                          
%    -----------| OBPF  |/                        
%         post  |       |\         
%         fiber  -------  \      Mach Zehnder     __|__
%                          \      ---------        / \
%                           \    /  +pi/4   \   /  ---     -------
%               		     \  /            \ /    |     |       | IRIC(:,2)                            
%                              -              -     |-----|  LPF  |-------
%                               \   ------   / \  __|__   |       |   
%                                \_| 1 bit|_/   \  / \     -------
%                                  | delay|        ---                 
%                                   ------          |
%
%   VARARGOUT has a variable number of arguments, from 1 to 2. In the
%   complete case it is VARARGOUT=[IRIC,SN], where Iric is a matrix
%   containing the received currents of channel ICH, while SN is a vector 
%   of the FFT coefficients of the electric field after the optical filter.
%   SN is used by BER_KL. 
%
%   X is a structure of fields:
%
%   X.oftype = optical filter (OBPF) type (see MYFILTER)
%   X.obw = OBPF 3 dB bandwidth normalized to the symbol rate. 
%   X.oord = optical filter order (for special filters, see MYFILTER)
%   X.eftype = electrical filter (LPF) type (see MYFILTER)
%   X.ebw = LPF 3-dB bandwidth normalized to the symbol rate. 
%   X.eord = electrical filter order (for special filters, see MYFILTER)
%   X.comp = component on which evaluate eye and calculate BER. Can be 
%            'phase' or 'quadrature'.
%
%   Optional parameters of X:
%
%   X.dpost = post compensating fiber cumulated dispersion [ps/nm]
%   X.slopez = post compensating fiber cumulated slope [ps/nm^2]
%   X.lambda = wavelength [nm] at which the post compensating fiber
%              has a cumulated dispersion equal to X.dpost.
%   X.b2b = 'b2b' evaluates the current in back-to-back configuration, i.e.
%              with the transmitter connected directly to the receiver. 
%              With this option the values of x.dpost and x.slopez are 
%              discarded.
%   X.mzdel = specify the delay of the upper brace of both the MZDI
%           interferometers of the DQPSK receiver. The default delay is 1 
%	    and it must be comprised in the interval 0 < mzdel <=1. Setting 
%	    mzdel to a value smaller than 1 implements the Partial DQPSK 
%	    (P-DQPSK) [1]. 
%
%   The post-compensating fiber is assumed as a purely ideal-linear fiber,
%   while the photodiodes are ideal (abs(.)^2).
%
%   Note: This function works over a copy of the electric field. All fields 
%   of the global variable GSTATE are left unchanged.
%
%   See also RECEIVER_OOK
%
%   [1] V. Mikhailov et al., "Experimental Investigation of Partial 
%   Demodulation of 85.3 Gb/s DQPSK signals," in. Proc. ECOC 2008, Bruxelles,
%   Belgium, 2008, paper We.1.E.5
%
%   Author: Paolo Serena, 2009
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

global CONSTANTS;  % CONSTANTS is a global structure variable.
CLIGHT = CONSTANTS.CLIGHT;      % speed of light in vacuum [m/s]

global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

if ~isfield(x,'oord')
    x.oord = 0;     % not using special filter
end

if ~isfield(x,'eord')
    x.eord = 0;     % not using special filter
end

if ~isfield(x,'mzdel')
    x.mzdel = 1;     % not using special filter
elseif (x.mzdel <= 0) || (x.mzdel > 1)
    error('error in ber_kl: mzdel must be  0 < mzdel <= 1')
end

%%%%%%%%%%% INITIALIZATION

Nfft = length(GSTATE.FN);
[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
npoints = 1:Nfft;                   % for all channels.
maxl=max(GSTATE.LAMBDA);
minl=min(GSTATE.LAMBDA);
lamc = 2*maxl*minl/(maxl+minl); %central wavelength: 1/lamc = 0.5(1/maxl+1/minl)
if nfc ~= GSTATE.NCH
    minfreq = GSTATE.FN(2)-GSTATE.FN(1);  
    deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich)); % frequency spacing    
    ndfn = round(deltafn./GSTATE.SYMBOLRATE/minfreq);  % spacing in points
    nind = nmod(npoints-ndfn,Nfft);    % left-circular shift    
    nch = 1;
    if ich == 1 % ndfnl & ndfnr are used for evaluating the energy after
        ndfnl = Nfft/2;
    else
        deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich-1)); % left channel
        ndfnl = round(deltafn./GSTATE.SYMBOLRATE/minfreq);
        ndfnl = round((ndfn-ndfnl)*0.5);    % spacing from left chan
    end
    if ich == GSTATE.NCH
        ndfnr = Nfft/2;
    else
        deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich+1)); % left channel
        ndfnr = round(deltafn./GSTATE.SYMBOLRATE/minfreq);        
        ndfnr = round((ndfnr-ndfn)*0.5);    % spacing from right chan
    end       
else
    nind = npoints;
    nch = ich;
    ndfnl = Nfft/2;
    ndfnr = Nfft/2;
end    
  
if isfield(x,'b2b')
    if strcmp(x.b2b,'b2b')
        b2b = 1;
        if isfield(x,'dpost'), x=rmfield(x,'dpost'); end;
    else
        error('the b2b field must be ''b2b''');
    end
else
    b2b = 0;
end

if b2b
    x.sigx = GSTATE.FIELDX_TX(:,nch);
else
    x.sigx = GSTATE.FIELDX(:,nch);
end

if isfield(x,'dpost')
    b20z = -x.lambda^2/2/pi/CLIGHT*x.dpost*1e-3; % beta2 [ns^2] @ lambda
    b30z = (x.lambda/2/pi/CLIGHT)^2*(2*x.lambda*x.dpost+...
        x.lambda^2*x.slopez)*1e-3;  
                                     % beta3 [ns^3] @ lambda

    % Domega_ik: [1/ns]. "i" -> at ch. i, "0" -> at lambda
    Domega_i0 = 2*pi*CLIGHT*(1./GSTATE.LAMBDA(ich)-1/x.lambda);    
    Domega_ic = 2*pi*CLIGHT*(1./GSTATE.LAMBDA(ich)-1/lamc);  
    Domega_c0 = 2*pi*CLIGHT*(1./lamc-1/x.lambda); 
    beta1z = b20z*Domega_ic+0.5*b30z*(Domega_i0^2-Domega_c0^2);    %[ns]    
    beta2z = b20z+b30z*Domega_i0;  % beta2*z [ns^2]@ GSTATE.LAMBDA
    % dispersion of the channels
    omega = 2*pi*GSTATE.SYMBOLRATE*GSTATE.FN';     % angular frequency [rad/ns]
    betat = omega*beta1z+0.5*omega.^2*beta2z+omega.^3*b30z/6;
    x.post_delay = GSTATE.SYMBOLRATE.*beta1z;

    Hf = fastexp(-betat);
else
    Hf = ones(Nfft,1);
    x.post_delay = 0;
end
Hf = Hf .* myfilter(x.oftype,GSTATE.FN,0.5*x.obw,x.oord);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% RECEIVER SIDE

x.sigx = fft(x.sigx);
x.sigx = x.sigx(nind);

x.avgebx = sum(abs(x.sigx(1:ndfnl)).^2) + sum(abs(x.sigx(Nfft-ndfnr+1:Nfft)).^2);
x.avgebx = x.avgebx/GSTATE.POWER(ich)/Nfft^2; % normalized average energy x bit
x.sigx = x.sigx .* Hf;                        % 1/3: optical filter

isy = ~isempty(GSTATE.FIELDY);
ndel = nmod((1:Nfft)-round(x.mzdel*GSTATE.NT),Nfft);           % 1 bit delay
if isy
    if b2b
        if  isempty(GSTATE.FIELDY_TX)
            x.sigy = zeros(nfr,1);
        else
            x.sigy = fft(GSTATE.FIELDY_TX(:,nch));
        end
    else
        x.sigy = fft(GSTATE.FIELDY(:,nch));
    end
    x.sigy = x.sigy(nind);
    x.avgeby = sum(abs(x.sigy(1:ndfnl)).^2) + sum(abs(x.sigy(Nfft-ndfnr+1:Nfft)).^2);
    x.avgeby = x.avgeby/GSTATE.POWER(ich)/Nfft^2; % normalized average energy x bit    
    x.sigy = x.sigy .* Hf;
    if nargout == 2    
        varargout(2) = {x};
    end
    x.sigx = ifft(x.sigx);
    x.sigy = ifft(x.sigy);
    
    Iric(:,1) = real(fastexp(-pi/4)*x.sigx.*conj(x.sigx(ndel)))+real(fastexp(-pi/4)*x.sigy.*conj(x.sigy(ndel)));     % 2/3: photodiode
    Iric(:,2) = real(fastexp(pi/4)*x.sigx.*conj(x.sigx(ndel)))+real(fastexp(pi/4)*x.sigy.*conj(x.sigy(ndel)));     % 2/3: photodiode   
else
    if nargout == 2    
        varargout(2) = {x};
    end
    x.sigx = ifft(x.sigx);

    Iric(:,1) = real(fastexp(-pi/4)*x.sigx.*conj(x.sigx(ndel)));     % 2/3: photodiode
    Iric(:,2) = real(fastexp(+pi/4)*x.sigx.*conj(x.sigx(ndel)));     % 2/3: photodiode
end
                                                % 3/3: lowpass filter
Hf = myfilter(x.eftype,GSTATE.FN,x.ebw,x.eord); % while the 3 dB bandwidth of the OBPF
Hf=repmat(Hf,1,2);
% goes from -x.obw/2 -> +x.obw/2, the LPF goes from 0 -> x.eftype.

    Iric = real(ifft(fft(Iric) .* Hf));
varargout(1) = {Iric};
