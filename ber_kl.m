function varargout=ber_kl(ich,x,pat)

%BER_KL Evaluate the ber for noncoherent transmission by Karhunen-Loeve method.
%   VARARGOUT=BER_KL(ICH,X,PAT) evaluates the bit error rate (ber) of
%   channel ICH of a non-coherent transmission by means of the Karhunen-Loeve 
%   (kl) method. Available modulation formats are OOK, DPSK, PSBT, DQPSK 
%   (see ELECTRICSOURCE).
%
%   PAT is the pattern for the ICH channel, after decoding (see 
%   PAT_DECODER).
%
%	X is a structure whose fields are:
%
%   X.rec = 'ook' to use receiver_ook, 'dpsk' to use receiver_dpsk,
%           'dqpsk' to use receiver_dqpsk  
%   X.oftype = optical filter type (see MYFILTER)
%   X.obw = optical filter 3 dB bandwidth normalized to the symbol rate
%   X.oord = optical filter order (for special filter, see MYFILTER)
%   X.eftype = electrical filter type (see MYFILTER)
%   X.ebw = electrical filter 3-dB bandwidth normalized to the symbol rate
%   X.eord = electrical filter order (for special filter, see MYFILTER)
%   X.osnr = Optical signal to noise ratios (osnr), [dB], over which the 
%            ber is evaluated. The osnr is over a conventional bandwidth 
%            of 0.1 nm and is measured immediately before the receiver. 
%            X.osnr refers to X.poln noise polarizations.
%   X.poln = Noise polarizations, 1 or 2. Note: X.poln is independent from 
%            the signal polarizations, e.g. the algorithm can work with two
%            noise polarizations and just one signal polarization. However,
%            with two signal polarizations X.poln must be 2.
%
%   X has also the following parameters required by the kl-method:
%
%   X.eta = bandwidth expansion factor. The kl method samples the signal
%           and the noise up to a frequency equal to X.eta times the  
%           bandwidth of the optical filter. Usually it is 1 < X.eta < 3. 
%   X.mu =  Time expansion factor. The memory of the receiver is X.mu times
%           the time duration of the memory devices inside the receiver,
%           i.e. the optical/electrical filter. For DPSK there is an 
%           additional memory due to the Mach-Zehnder delay interferometer.
%           The memory of such devices is approximated by the inverse of 
%           their bandwidths, as suggested in [1]. 
%           Usually it is 1 < X.mu < 10. 
%   X.saddle = 'yes': the ber is evaluated through the saddle point
%           approximation (faster). 'no': the ber is evaluated by numerical
%           integration of the moment generating function (slower, but more 
%           accurate).
%
%   For more details about X.eta, X.mu and X.saddle see [1].
%   X can also have the optional parameters:
%
%   X.ber = reference ber at which the algorithm returns the corresponding
%           osnr, searched within the range X.osnr by numerical 
%           interpolation.
%   X.interp = interpolation method for finding X.ber, see INTERP1.
%           Default is 'spline'.
%   X.extrap = 'yes': X.ber can be extrapolated outside X.osnr, see
%           INTERP1. 'no': If X.ber is outside the range X.osnr the 
%           function returns VARARGOUT(2) = NaN, which is also the default 
%           strategy.
%   X.plot = 'ploteye': plots the eye in the active figure; 'plotcur' plots
%           the received current.
%   X.color = color string for the plot (see PLOT). E.g. 'b-'.
%   X.dpost = post compensating fiber cumulated dispersion [ps/nm]
%   X.slopez = post compensating fiber cumulated slope  [ps/nm^2]
%   X.lambda = wavelength [nm] at which the post compensating fiber
%           has a cumulated dispersion equal to X.dpost.
%   X.comp = component on which evaluate eye and calculate BER (DQPSK
%            modulation only). Can be 'phase' or 'quadrature'.
%   X.b2b = 'b2b': The function works in back-to-back transmission.
%   X.print = structure for print. E.g. X.print = {'nomefile','eye'} or
%            X.print = {'nomefile','current'}, prints to file nomefile the 
%            eye or the current, respectively. nomefile will be place into
%            GSTATE.DIR within a directory ending with '.MOD'.
%   X.delay = 'theory' means that the delay uses the theoretical delay saved 
%           within GSTATE.DELAY (see CREATE_FIELD) Per default the delay is 
%           measured by a cross-correlation measurement between the 
%           received current and an artificial pulse amplitude modulation 
%           (PAM) signal with ideal non-return to zero bits with symbols 
%           equal to PAT. The correlation method is useful in presence of
%           polarization mode dispersion (PMD).
%	X.threshold = Fixed threshold for the threshold detector. By default the 
%			threshold is optimized for OOK and set to zero for DPSK/DQPSK.
%
%   X.mzdel = specify the delay of the upper brace of the MZDI
%           interferometer for DPSK/DQPSK. The default delay is 1 and it
%           must be comprised in the interval 0 < mzdel <=1. Setting mzdel
%           to a value smaller than 1 implements the Partial DPSK/DQPSK 
%           (P-DPSK/DQPSK) [3]. 
%
%
%   The ook receiver is composed of an ideal, purely linear, post compensating 
%   fiber + optical filter + photodiode + electrical lowpass filter 
%   (see RECEIVER_OOK or RECEIVER_DPSK). For DPSK there is also a
%   Mach-Zehnder interferometer before the photodiodes.
%
%   VARARGOUT can contain a variable number of argument, from 1 to 3. In the
%   complete case it is VARARGOUT=[PB,OSNR,EO] where PB is a vector of
%   size(X.osnr) containing the ber corresponding to X.osnr, OSNR is
%   the signal-to-noise ratio [dB] that yields X.ber while EO is the 
%   normalized eye opening (see EVAL_EYE).
%
%   Note 1: This function works over a copy of the electric field. All  
%           fields of the global variable GSTATE are left unchanged.
%   Note 2: all fields of X must be lowercase.
%   Note 3: The noise is assumed white over the frequency. The possible
%       presence of parametric gain is not accounted by this function.
%
%   See also PATTERN, MYFILTER, RECEIVER_OOK, RECEIVER_DPSK, BEST_EYE
%            EVAL_EYE, BEST_SP, PAT_DECODER
%
%   This function implements the algorithm proposed by E. Forestieri in [1]
%
%   [1] E. Forestieri, "Evaluating the Error Probability in Lightwave 
%       Systems with Chromatic Dispersion, Arbitrary Pulse Shape and Pre-  
%       and Postdetection Filtering,", J. Lightw. Technol., vol. 18, no.11,
%       pag. 1493-1503, Nov. 2000.
%
%   The DPSK version of the algorithm can be found in [2]. Note that in
%   this function noise parametric gain is neglected.
%
%   [2] P. Serena, A. Orlandini and A. Bononi, "A Parametric-Gain Approach
%       to the Analysis of Single-Channel DPSK/DQPSK Systems With Nonlinear
%       Phase Noise," J. Lightw. Technol., vol. 24, no. 5, pag. 2026-2037,
%       May 2006.
%
%   [3] B. Mikkelsen, C. Rasmussen, P. Mamyshev and F. Liu, "Partial DPSK 
%       with excellent ﬁlter tolerance and OSNR sensitivity," Electronics
%       Letters, vol. 42, no. 23, pp. 1363-1364 , Nov. 2006
%
%   Many thanks to E. Forestieri for thefortran code of his algorithm, 
%   which has been the main inspiration for this function.
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
global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

%%%%%%%%%%%%%%%
CLIGHT = CONSTANTS.CLIGHT;  % speed of light [m/s]
INTMETH = 'spline'; % type of interpolation for finding OSNR @ x.ber
EXTRVAL = NaN; % default: do not extrapolate OSNR @ x.ber outside the range
%%%%%%%%%%%%%%%

%%%%% Initialize and check

[Nfft,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
fnames = fieldnames(x);
isy = ~isempty(GSTATE.FIELDY);

if any(strcmp(fnames,'saddle')) % same as isfield (faster)
    if strcmp(x.saddle(1),'y')
        sad = 1;
    else
        sad = 0;
    end
else
    sad = 0;
end

if any(strcmp(fnames,'interp')) % same as isfield (faster)
    meth = x.interp;
else
    meth = INTMETH;
end
if any(strcmp(fnames,'extrap')) % same as isfield (faster)
    if strcmp(x.extrap(1),'y')
        extr = 'extrap';
    else
        extr = x.extr;
    end
else
    extr = EXTRVAL;
end
if nargout == 0
    error('missing output arguments');
end
if ~any(strcmp(fnames,'oord')) % same as isfield (faster)
    x.oord = 0;     % not using special filter
end

if ~any(strcmp(fnames,'eord')) % same as isfield (faster)
    x.eord = 0;     % not using special filter
end

if ~isfield(x,'mzdel')
    x.mzdel = 1;     % not using special filter
elseif (x.mzdel <= 0) || (x.mzdel > 1)
    error('error in ber_kl: mzdel must be  0 < mzdel <= 1')
end

%%%%% Set the bandwidth

if strcmp(x.oftype,'movavg')
    Bo = x.obw;     % here 1/x.obw is the duration of 
                    % the mov. average.
    Bon = x.obw;    % Because with ideal OOK, x.obw=1 (matched filter), 
    % and no LPF, I want just one noise sample, i.e., M=1 and T0=1.
    x.obw=Bo*2;     % Because he receiver thinks that Bo is a bandpass bandwidth
else
    Bo = x.obw/2;   % lowpass bandwidth
    Bon = b3dB2neq(x.oftype,Bo,x.oord);  % noise equivalent bandwidth
end
Brn = b3dB2neq(x.eftype,x.ebw,x.eord); % noise eq. bandwidth of electrical filter
if strcmp(x.eftype,'movavg')
    if x.mu > 1
        Brn = x.ebw*x.mu; % the duration is exactly 1/x.ebw. Mult. by x.mu 
    end                   % in order to save an if next.
end
pat = pat(:);
%%%%% Set the number of frequencies

if x.eta > 0
    LE = ceil(x.eta*GSTATE.NSYMB*Bon);
%    LL = min(Nfft/2,2^nextpow2(LE))-1; % there are only Nfft frequencies available!
    dLE = ceil(2*LE/GSTATE.NSYMB)*GSTATE.NSYMB; % make 2*LE a multiple of GSTATE.NSYMB
    LE = ceil(dLE/2); % make LE an even number (necessary for matrix V, see next).
    if LE > Nfft/2
        LE = Nfft/2;
        x.eta = LE/(GSTATE.NSYMB*Bon);
        warning('optilux:berkl',['x.eta too large. Reduced to ',num2str(x.eta)]);
    end    
else
    LE = 0; % get only the DC component
end

%%%%% Get the electric field FFT (sn),the sampled current (Iric), and the
%%%%% sampling time (tk).

[eo,ts,Iric,sn,avgdelay] = eval_eye(ich,x,pat); % get the current
Iric = Iric/GSTATE.POWER(ich);   % normalize to transmitted peak power

if strcmp(x.rec,'dpsk') || (strcmp(x.rec,'dqpsk') && strcmp(x.comp,'phase'))
    sgn = 1-2*pat;
    if any(strcmp(fnames,'threshold'))
        gthend = x.threshold; % user defined threshold
        gthlow = x.threshold;
    else
        gthend = 0;
        gthlow = 0;
    end
elseif (strcmp(x.rec,'dqpsk') && strcmp(x.comp,'quadrature'))
    sgn = 1-2*pat;
    if any(strcmp(fnames,'threshold'))
        gthend = x.threshold; % user defined threshold
        gthlow = x.threshold;
    else
        gthend = 0;
        gthlow = 0;
    end
elseif strcmp(x.rec,'ook')
    sgn = 1-2*pat;
    if any(strcmp(fnames,'threshold'))
        gthend = x.threshold; % user defined threshold
        gthlow = x.threshold;
    else
        gthend = max(Iric);
        gthlow = 0;
    end
end

%%%%% Delay estimation


tk = ts + avgdelay;   % sampling time (avgdelay is a double)

%%%%% Set the memory of the KL method

if strcmp(x.oftype,'ideal')
    Bopeta = Bon;
else
    Bopeta = x.eta*Bon;
end
if strcmp(x.oftype,'movavg') && (x.mu > 1)
    T0 = 1/Bon + x.mu/Brn;    % memory time
else
    T0 = x.mu*(1/Bon + 1/Brn);
end
if strcmp(x.rec,'dpsk') || strcmp(x.rec,'dqpsk'), T0 = T0 + 1;end; % +1: memory of MZ

M = ceil(Bopeta*T0);   % The number of noise samples is 2*M+1
M2p1 = 2*M+1;

%%%%% Create matrix A, its eigenvectors and eigenvalues.

fnos = 0:1/T0:2*M/T0;  % freq. of (A.12) in [1]
fnts = fnos - M/T0;

qfirstrow = myfilter(x.eftype,fnos,x.ebw,x.eord).'; 
Qmat = toeplitz(qfirstrow);             % (A.12) 
Ho = abs(myfilter(x.oftype,fnts,Bo,x.oord));
Hmat = diag(Ho);   % (A.17)
if strcmp(x.rec,'dpsk')
    dang = -2*pi*(-M:M)/T0*x.mzdel;
    dddiag = complex(cos(dang),sin(dang));
    Dmul = diag(dddiag); % see point 3) in Section V of [2], no pg case.
    Qmat = 0.5*(Dmul'*Qmat + Qmat'*Dmul);
elseif strcmp(x.rec,'dqpsk')
    dang = -2*pi*(-M:M)/T0*x.mzdel;
    dddiag = complex(cos(dang),sin(dang));
    Dmul = diag(dddiag); % see point 3) in Section V of [2], no pg case.
    if strcmp(x.comp,'phase')
        Qmat = 0.5*(Dmul'*Qmat*fastexp(-0.25*pi) + Qmat'*Dmul*fastexp(0.25*pi));
    else
        Qmat = 0.5*(Dmul'*Qmat*fastexp(0.25*pi) + Qmat'*Dmul*fastexp(-0.25*pi));
    end
else
    dddiag = []; % not used
end
Amat = Hmat'*Qmat*Hmat; % (A.20)
Amat = 0.5*(Amat+Amat');  % remove numerical errors of order 10^-15.

[Umat,Dmat] = eig(Amat);
ld = real(diag(Dmat));    % eigenvalues
ldex(1) = min(ld);
ldex(2) = max(ld);
ld2 = ld.^2;
if find(ld == 0)
    error('Cannot evaluate the BER: one or more zero eigenvalues.');
end

%%%%% Create matrix V

V = create_V(ich,tk,Nfft,LE,sn.sigx,M2p1,M,T0,x,dddiag);
Eb = sn.avgebx; % normalized average energy x bit
if isy
    if x.poln == 1
        error(['two signal polarizations but just one ',...
            'noise polarization']);
    end
    Vy = create_V(ich,tk,Nfft,LE,sn.sigy,M2p1,M,T0,x,dddiag);
    Eb = Eb + sn.avgeby;
end

%%%%% Set up MGF parameters

b2 = zeros(M2p1,GSTATE.NSYMB);
for n=1:M2p1
    b2(n,:) = Ho(n)*V(:,n).';  % H^(T*)*V in (A.26) of [1]
end
b2 = Umat' * b2;   % (A.26) of [1]. Note: here b2 is a matrix.
b2 = abs(b2).^2;	% |b_i|^2. Matrix (M2p1 x Nsymb) 
if isy
    b2y = zeros(M2p1,GSTATE.NSYMB);
    for n=1:M2p1
        b2y(n,:) = Ho(n)*Vy(:,n).';  % H^(T*)*V in (A.26) of [1]
    end
    b2y = Umat' * b2y;      % (A.26) of [1]. Note: here b2 is a matrix.
    b2 = b2 + abs(b2y).^2;	% |b_i|^2. Matrix (M2p1 x Nsymb)
end
sl = sum(ld);
sl2 = sum(ld2);
sb2 = sum(b2);
uld = 1./ld;
sb2l = uld' * b2;  % (15) of [1]. Vector (1 x Nsymb)
sb2l = sb2l(:);

%%%%% Evaluate the BER

cOSNR = CLIGHT/GSTATE.LAMBDA(ich)^2/GSTATE.SYMBOLRATE*0.1; % OSNR is over 0.1 nm!
LOSNR = length(x.osnr);
pb = ones(1,LOSNR);         % error probability
for nosnr=1:LOSNR
    rho = 10^(x.osnr(nosnr)*0.1);
    sg2 = Eb/(2*x.poln*cOSNR*rho*T0);   % noise variance
    dsg2 = 2*sg2;
    qsg2 = 4*sg2;
    qsg4 = dsg2*dsg2;
    
    % Now search the optimal threshold using fminbnd.
    [gthopt,pb(nosnr)] = fminbnd(@eval_ber,gthlow,gthend,[],Iric,b2,sb2,...
        sb2l,sgn,sg2,dsg2,qsg2,qsg4,sl,sl2,x.poln,ld,ldex,ld2,M2p1,sad);
    
end % for nosnr=1:LOSNR

%%%%%% Output

varargout(1)={pb};
if nargout >= 2    
    indf = isfinite(pb);
    rngosnr = x.osnr(indf);
    if length(rngosnr) >= 2
        pblog = log10(pb(indf)); % the interpolation is in a log10 scale
        pblog = pblog(pblog ~= 0);
        if any(strcmp(fnames,'ber')) % same as isfield (faster)
            osnrrif = interp1(pblog(isfinite(pblog)),rngosnr(isfinite(pblog)),log10(x.ber),meth,extr);
        else
            osnrrif = NaN;
        end
        varargout(2) = {osnrrif};
    else
        varargout(2) = {NaN};
        warning('optilux:ber_kl','There should be at least two OSNR.');
    end
end
if nargout == 3
    varargout(3)={eo};
end

if GSTATE.PRINT

    %%%% PRINT summary
    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===              ber_kl              ===\n');
    fprintf(fid,'========================================\n\n');
	if any(strcmp(fnames,'b2b')), fprintf(fid,'BER in back-to-back\n');end
	fprintf(fid,'\neta factor: %.2f. mu factor: %.2f. Memory: %.3f\n',...
		x.eta,x.mu,T0);
	fprintf(fid,'noise polarizations: %d\n',x.poln);
	if sad == 0
  		fprintf(fid,'BER measured with numerical integration \n');
	else
  		fprintf(fid,'BER measured with saddlepoint approximation\n');
	end	
	fprintf(fid,'Receiver type: %s\n',x.rec);
	fprintf(fid,'Optical filter type: %s\n',x.oftype);
	fprintf(fid,'Optical filter bandwidth: %.2f\n',x.obw);
	if x.oord ~= 0
  		fprintf(fid,'Optical filter order: %d\n',x.oord);
	end
	fprintf(fid,'Electrical filter type: %s\n',x.eftype);
	fprintf(fid,'Electrical filter bandwidth: %.2f\n',x.ebw);
	if x.eord ~= 0
  		fprintf(fid,'Electrical filter order: %d\n',x.eord);
	end
	if any(strcmp(fnames,'dpost'))
		fprintf(fid,'\nPost compensation: %.3f [ps/nm] @lambda=%.2f [nm]\n',...
		x.dpost,x.lambda);
	end    
    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function V = create_V(ich,tk,Nfft,LE,sn,M2p1,M,T0,x,dddiag)

%%%%% Create matrix V
%
% The columns of matrix V are the vector v, (A.13) of [1]. A fast way for
% implementing such matrix is described in [2]. Strictly speaking, the
% matrix V satisfies V=F*Ih*G, where Ih=[I,I,...,I] is the GSTATE.NT
% repetition of I=eye(GSTATE.NSYMB), F is the Fourier matrix, G the matrix 
% \hat{S} in [2], Section 5. The first product can be done by the ifft
% algorithm, the second involves only GSTATE.NT non-zero multiplications.

global GSTATE

dtn = 2i*pi*tk;
fred = [GSTATE.FN(1:LE),GSTATE.FN(Nfft-LE+1:Nfft)]; % freq up to LE
% index l in (A.1) of [1] is -LE <= l <= LE-1  
ldown = 1:GSTATE.NSYMB:2*LE;  % down-sampling
% Note: snred is the signal over fred normalized to peak Tx power
snred = [sn(1:LE); sn(Nfft-LE+1:Nfft)]/Nfft/sqrt(GSTATE.POWER(ich));   
V = zeros(GSTATE.NSYMB,M2p1);

for n=1:M2p1
    fkn = fred - (n-M-1)/T0;
    allHf = myfilter(x.eftype,fkn,x.ebw,x.eord).'; % row vector        
    for k=0:GSTATE.NSYMB-1
        lind = k+ldown;               % convert to fftshift notation
        Hf = exp(dtn*fred(lind)).*allHf(lind);
        V(k+1,n) = Hf*snred(lind);    % scalar product
    end
end
V = ifft(V)*GSTATE.NSYMB; % multiplication by the Fourier matrix -> ifft

if strcmp(x.rec,'dpsk') % see point 4) in Section V of [2], no pg case.
    dddiag = conj(dddiag);
    Dmul_dagger = repmat(dddiag,GSTATE.NSYMB,1);
    V = 0.5*(Dmul_dagger.*V+circshift(V,1));  
elseif strcmp(x.rec,'dqpsk') % see point 4) in Section V of [2], no pg case.
    dddiag = conj(dddiag);
    Dmul_dagger = repmat(dddiag,GSTATE.NSYMB,1);
    if strcmp(x.comp,'phase')
        V = 0.5*(Dmul_dagger.*V*fastexp(-0.25*pi)+circshift(V,1)*fastexp(0.25*pi));
    else
        V = 0.5*(Dmul_dagger.*V*fastexp(0.25*pi)+circshift(V,1)*fastexp(-0.25*pi));        
    end
end

%--------------------------------------------------------------------------
function pb=eval_ber(gth,Iric,b2,sb2,sb2l,sgn,sg2,dsg2,qsg2,qsg4,...
    sl,sl2,pol2,ld,ldex,ld2,M2p1,sad)

% Returns the BER pb given the threshold gth.

global GSTATE

xi = gth - Iric + sb2l; % (22) of [1]

%%%%% Find the saddlepoint 

eta = dsg2*sl*pol2+sb2l;           % (19) of [1]
vars = qsg2*(sg2*sl2*pol2+sb2);    % (20) of [1]
%            var_tot(k) = vars; % sampled noise variance
%            var_ssp(k) = qsg2*sb2(k);  % signal-noise variance
dx = xi-eta;
% u0 is the saddlepoint for each bit: vector [GSTATE.NSYMB,1]
[u0,nsw] = saddle(sgn,xi,dx,dsg2,qsg2,qsg4,vars,ld,ldex,...
    ld2,b2,pol2);  % roots of (B.8) of [1]

if any(nsw)
    fnsw1 = find(nsw == 1);
    fnsw2 = find(nsw == 2);
    if ~isempty(fnsw1)
        warning('optilux:ber_kl',[' no saddle point at bit: ',...
            num2str(fnsw1)]);
    end
    if ~isempty(fnsw2) 
        warning('optilux:ber_kl',[' inaccurate saddle point at bit: ',...
            num2str(fnsw2)]);
    end
end
%%%%% Evaluate the derivatives of \Phi_{n_k}

u02 = u0.*u0;
ds2u0 = dsg2*u0;
ds2u02 = dsg2*u02;
bbs = ld*ds2u0; % matrix whose kth column is ld*ds2u0(k)
umbs = 1-bbs;
dab2 = qsg2*b2;
repld2 = repmat(ld2,1,GSTATE.NSYMB);
b2umbs = qsg4*repld2.*umbs;
dn2 = dab2+b2umbs*pol2;
dd2 = umbs.*umbs.*umbs;
fi2 = 1./u02 + sum(dn2./dd2);        % (B.10) of [1]
repds2u02 = repmat(ds2u02,M2p1,1);
psr = sum(repds2u02.*b2./umbs-pol2*(log(umbs)+bbs)); % (B.18) of [1]
                                                     % for all bits.

% ek is a vector containing the error probability of the bits.
if sad                  % saddlepoint approximation
    ek = sgn.'.*exp(-dx.'.*u0+psr)./(u0.*sqrt(2*pi*fi2));
else                    % numerical inversion of the MGF
    repld = repmat(ld,1,GSTATE.NSYMB);
    dn3 = dsg2*repld.*(3*dab2+2*b2umbs*pol2);
    dd3 = dd2.*umbs;
    fi3 = -2./(u02.*u0) + sum(dn3./dd3);  % (B.14) of [1]
    dw = sqrt(0.5./fi2); % (B.17) of [1]. Step of integration
    cur = fi3./(3*fi2);  % (B.13) of [1]. Parabola curvature
    cur(cur < 0) = 0;
    hc = 0.5*cur;

    [ek,lowacc] = intmgf(sgn,dx,u0,cur,psr,dw,hc,dsg2,ld,...
        b2,pol2);
    
    if any(lowacc)
        warning('optilux:ber_kl',['inaccurate BER during the',...
            ' search of the best gth @ gth=',num2str(gth)]); 
    end
end

pb = sum(ek);
pb = pb/GSTATE.NSYMB;    % error probability
pb(imag(pb) ~= 0 || pb < 0 || pb > 1) = NaN;

%--------------------------------------------------------------------------    
function [ek,lowacc] = intmgf(sgn,dx,u0,cur,psr,dw,hc,dsg2,ld,b2,pol2)
 
% Evaluate the BER by numerical integration of the MGF. See [1] eq. (B.15)

%%%%%%%%%%%%%%%
TOLEK = 1e-4;   % tolerance of the integral.
TOLEPS = 1e-6;  % tolerance of the function, see comment after (B.16) in [1]
MAXEPS = 100;   % max value of eps2 for using the parabola (B.12) in [1]
MAXITER = 10;   % maximum number of recursive iterations.
%%%%%%%%%%%%%%%

nsymb = length(sgn);
lowacc = zeros(1,nsymb); % 1: low accuracy, 0: ok.
ek = ones(1,nsymb);
for k=1:nsymb
    
    ek(k) = 0.5*sgn(k)*exp(-dx(k)*u0(k)+psr(k))/u0(k);    % f(0)/2 in (B.16) of [1]
    ekp = ek(k);
    w0 = 0;     % Shift. It is w0 ~= 0 when the recursion is on, and the 
                % trapezoidal method adds point on a shifted grid. 
    w1 = dw(k);    % Step of integration         
    nv = 1;     % Number of trapezoidal recursions. If greater than MAXITER,  
                % the result may be of low accuracy.
    nx = 0;

    % evaluate the integral (B.16) in [1] by trapezoidal rule. The integration
    % is recursive by halving the step until the relative error is reached. 
    cond = 1;
    while cond

        % The number of function evaluation in the integral is stopped when
        % f(n\Delta v) in [1], eq. (B.16)  is negligible, i.e. eps < TOLEPS.

        eps2 = Inf;  % Relative error on the function.
        while (eps2 > TOLEPS)                                    
            nx = nx+1;
            dwn = w0+nx*dw(k); % (n\Delta v) in (B.16) of [1]
            un = u0(k)+hc(k)*dwn*dwn + i*dwn; % (B.12) of [1]: parabola
            % un is a good approx. of the steepest descent path
            ds2un = dsg2*un;
            ds2un2 = ds2un*un;
            bis = ds2un*ld;
            umbis = 1.d0-bis;
            psif = sum(ds2un2*b2(:,k)./umbis-pol2*(log(umbis)+bis));
            psif = (1-i*cur(k)*dwn)*exp(-dx(k)*un+psif)/un; %f(n\Delta v)
            eps2 = abs(psif/(dw(k)*ek(k))); % relative error
            ek(k) = ek(k)+sgn(k)*real(psif);  % (B.16) of [1]
            if (eps2 > MAXEPS)  % repeat over a straight line path
                cur(k) = 0;
                hc(k) = 0;
                ek(k) = 0.5*sgn(k)*exp(-dx(k)*u0(k)+psr(k))/u0(k);    
                ekp = ek(k);
                w0 = 0;     
                nv = 1;  
                w1 = dw(k);
            end
        end
        cond = abs((ekp-0.5*ek(k))/ekp) > TOLEK;   % reached the tolerance?
        if (nv < MAXITER)
            ekp = ek(k);
            w0 = 0.5*w1;    % half the step
            dw(k) = w1;
            w1 = w0;
            nx = -1;
            nv = nv+1;
        else
            lowacc(k) = 1; % low accuracy. exit and print a warning
    %        warning('optilux:berkl','Warning: low accuracy in ber_kl.m');
            break
        end
    end
    if (nv > 1), dw(k) = 2*w0; end;
    ek(k) = ek(k)*dw(k)/pi; % error probability

end

% ------------------------------------------------------------------------

function beqn=b3dB2neq(ftype,bw,ord)

%B3DB2NEQ Convert the 3 dB bandwidth into the noise equivalent bandwidth.
%   BEQN=B3dB2NEQ(FTYPE,BW,ORD) returns in BEQN the noise equivalent 
%   bandwidth of the filter FTYPE having the 3dB bandwidth BW. FTYPE can be
%   one of the following string:
%
%   'movavg'     : Short term integrator (moving average) (see Note 1)
%   'gauss'      : Gaussian filter 
%   'butt2'      : Butterworth 2nd order
%   'butt4'      : Butterworth 4th order
%   'butt6'      : Butterworth 6th order
%   'ideal'      : Ideal filter
%   'bessel5'    : Bessel 5th order
%   'rc1'        : RC1 filter
%   'rc2'        : RC2 filter
%   'supergauss' : Super-Gaussian of order ORD
%
%   Note 1: For the moving average BW is not the 3dB bandwidth,
%           but the first zero of the sinc, i.e. 1/BW is the time-duration 
%           of the moving average. The 3dB bandwidth is 0.443*BW.
%
%   Thanks to E. Forestieri for providing the filter expressions.
%
%   Author: Paolo Serena, 2007
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2007  Paolo Serena, <serena@tlc.unipr.it>
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


% CONSTANTS

      r2m1=0.41421356237309504880;
      r2p1=2.41421356237309504880;


switch ftype
    
  case 'movavg'          
% Short term integrator

      beqn = 0.5*bw;
      
  case {'gauss','gauss_off'}
% Gaussian

      beqn = 0.5*sqrt(pi/log(2))*bw;

  case 'butt2'
% Butterworth 2nd order

      beqn = bw/sinc(0.25);
      
  case 'butt4'
% Butterworth 4th order

      beqn = bw/sinc(0.125);
      
  case 'butt6'
% Butterworth 6th order

      beqn = bw/sinc(1/12);
      
  case 'ideal'
% Ideal filter

      beqn = bw;
      
  case 'bessel5'
% Bessel 5th order

      beqn = 1.039*bw;
      
  case 'rc1'
% RC filter 1st order

      beqn = 0.5*pi*bw;

  case 'rc2'
% RC filter 2nd order

      beqn = 0.25*pi*bw/(r2m1*sqrt(r2p1));
      
  case 'supergauss'    
% Super-Gaussian of order ORD
      if ord == 0
          error('supergauss order must be > 0');
      end
      udo = 1/(2*ord);
      beqn = bw*gamma(1+udo)*log(2)^(-udo);

  otherwise
      error('the filter ftype does not exist.');      
end
