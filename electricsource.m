function elec = electricsource(pat, format, symbrate, ptype, duty, roll, varargin)
%ELECTRICSOURCE creates the electric modulating signal
%   ELEC = ELECTRICSOURCE(PAT, FORMAT, SYMBRATE, PTYPE, DUTY, ROLL, INPOW, ORD,PAR)
%   creates the electric signal which is passed to one input of the  
%   modulator, using the pattern PAT created by PATTERN. 
%   FORMAT is a string  that describes the modulation format. 
%
%   SYMBRATE is the signal's baudrate in [Gbaud], and it is associated to
%   the global variable GSTATE.SYMBOLRATE.
%
%   PTYPE is the pulse type and it is one of the following strings:
%
%   <string used in MYFILTER>:      Filters an ideal signal with the
%                                   correspondent filter in MYFILTER. In
%                                   this case ROLL is the 3dB bandwidth and
%                                   ORD is the order for special filters.
%
%   'cosroll':  Pulses with a raised cosine behavior during the
%               commutation states. In this case 0<ROLL<=1 indicates the
%               roll-off. The elementary pulse assumes the form [1]:
%                _
%               |    1       0 <= abs(t) <= (1-roll)/2*duty
%               |
%       p(t) = .     1/2*{1+cos[pi/roll/duty*...
%               |        (abs(t)-(1-roll)*duty/2)]},
%               !         if  (1-roll)*duty/2 <= abs(t) <= (1+roll)*duty/2
%               |
%                _   0      abs(t) > (1+roll)/2*duty
%
%   'dirac':      Dirac's delta pulses.
%   'idpulse':    Ideal pulses with only two levels. Do not confuse with
%                 the string 'ideal' which calls for the ideal filter in
%                 MYFILTER.
%   'sech':       Bright solitons (still to be implemented)
%   'tanh':       Dark solitons (still to be implemented)
%
%   ROLL must be always declared. If you don't need ROLL, set it, for
%   instance, equal to the empty variable, i.e. ROLL=[].
%
%   INPOW (optional) is 'power' if the shaping is on the signal's power
%   (abs(.)^2), otherwise the shaping is over the electric field (abs(.)).
%   DUTY is the duty-cycle, and must be 0 <= DUTY <= 1. 
%
%   Supported modulation formats (in square brackets the required modulator) 
%   include:
%       - OOK [MZ_MODULATOR]
%       - BPSK (and thus DPSK and NF-DPSK) [MZ_MODULATOR]
%       - PSBT [MZ_MODULATOR] In this case you can also specify two fields in the 
%         structure of parameters PAR to model the electrical filter:
%         * efilt    = electrical filter type (default: bessel5)
%         * efiltbw  = electrical filter bandwidth (default: 0.3)
%       - QPSK (and thus DQPSK and NF-DQPSK). Only the driving signal for
%         one of the two quadratures is created and thus ELECTRICSOURCE
%         must be called two times to get the required inputs of
%         QI_MODULATOR
%       - USERDEF In this case custom electric signal are generated: the
%         user must specify the following field of the structure PAR:
%         * alphabet = size of the alphabet of PAT
%         * limits   = a 2x1 vector containing lower and upper values of
%           the generated signal. Symbols are assumed equally spaced
%         OR
%         * ampls    = a ALPHABETx1 vector containing amplitudes associated
%           to every symbol of PAT
%
%   See also PATTERN, MZ_MODULATOR, LASERSOURCE
%
%   Author: Marco Bertolini, 2009
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

global GSTATE;

if length(symbrate) == 1
    GSTATE.SYMBOLRATE = symbrate;    % this is only a global initialization,
                                    % since the symbrate is not used in this function
else
    error('symbrate must have length 1');
end

if ~ischar(format)
    error('The modulation format must be a string')
end

for arg=1:length(varargin)
    if isnumeric(varargin{arg})
        ord = varargin{arg};
    elseif ischar(varargin{arg})
        inpow = varargin{arg};
    elseif isstruct(varargin{arg})
        par = varargin{arg};
    else
        error('inpow, ord or par wrong format')
    end
end
if ~exist('inpow','var')
    inpow='';
end
if ~exist('ord','var')
    ord=0;
end
if ~exist('par','var')
    par = [];
end

if ~ischar(ptype)
    error('ptype must be a string');
end


% prepare elctric fields
elec = zeros(GSTATE.NSYMB*GSTATE.NT,1);

if strcmp(format,'psbt')
    par=checkparam(par,format);
end

switch format
    case {'ook'}
        if max(pat) > 1
            error('%s requires a binary pattern',format);
        end
        elec = elecsrc(pat,ptype,duty,roll,inpow,ord);
    case {'bpsk','nf-dpsk','qpsk','nf-dqpsk'}
        if max(pat) > 1
            error('%s requires a binary pattern',format);
        end        
        elec = elecsrc(2*pat-1,ptype,duty,roll,inpow,ord);
        if strcmp(format,'nf-dpsk') || strcmp(format,'nf-dqpsk')
            elec=circshift(elec,-GSTATE.NT/2);
            % a shift of half symbol period is present in
            % all "narrow filtered" modulation formats (suffix -nf)
        end
    case {'psbt'}
        if max(pat) > 1
            error('%s requires a binary pattern',format);
        end        
        elec = elecsrc(2*pat-1,ptype,duty,roll,inpow,ord);
        elec = lpfilter(elec,par.efilt,par.efiltbw);        
        elec=circshift(elec,-GSTATE.NT/2);
    case 'userdef'
        if ~isfield(par,'alphabet') || ~xor(~isfield(par,'limits'),~isfield(par,'ampls'))
            error(strcat('when using ''userdef'' format you must provide ',...
                         'the alphabet and either the limits or the amplitude values'))
        end
        if isfield(par,'ampls')
            if par.alphabet ~= length(par.ampls)
                error('the number of amplitudes must be equal to the size of the alphabet');
            end            
            ampl = zeros(size(pat));
            for alph = 0:par.alphabet-1
                ampl(pat==alph)=par.ampls(alph+1);
            end            
        else
            if length(par.limits) ~= 2
                error('you must provide lower and upper limits for the amplitude');
            end            
            if par.limits(1) > par.limits(2)
                error('you must provide the amplitude limits in ascending order');
            end                        
            dim = par.alphabet;
            limits = par.limits;
            step = (limits(2)-limits(1))/(dim-1);
            ampl = limits(1)+step*pat;
        end
        elec = elecsrc(ampl,ptype,duty,roll,inpow,ord);
    otherwise
        error('wrong modulation format');
end


%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function elec = elecsrc(pat,ptype,duty,roll,inpow,ord)

global GSTATE;

if max(duty) > 1
    error('duty must be <= 1');
end
if exist('roll','var') ~= 1
    error('Undefined roll');
end
if strcmp(ptype,'cosroll') && (max(roll) > 1)
    error('roll-off must be <= 1');
end

if ~(strcmp(ptype,'cosroll') || strcmp(ptype,'dirac') || strcmp(ptype,'idpulse')...
        || strcmp(ptype,'sech') || strcmp(ptype,'tanh'))

    flag = 1;   % after modulation, apply the filter.
    ptypet = 'idpulse';     % temp ptype
else
    flag = 0;   % after modulation, exit the function.
    ptypet = ptype;
end


% Modulate

Nfft = GSTATE.NT*GSTATE.NSYMB;               % number of points
elec = zeros(Nfft,1);

elpulse = pulseshape(ptypet,roll,duty,GSTATE.NT);  % single pulse

nstart = Nfft-GSTATE.NT+1;      % The first elpulse starts at the end of
                                % the sequence (cyclic periodicity).
nend = Nfft;
elec(nstart:nend) = pat(1)*elpulse(1:GSTATE.NT);
elec(1:GSTATE.NT) = pat(1)*elpulse(GSTATE.NT+1:GSTATE.NT*2);
for kbit=2:GSTATE.NSYMB          % all other bits: 2 -> Nsymb
    nstart = (kbit-2)*GSTATE.NT+1;
    nend = kbit*GSTATE.NT;
    elec(nstart:nend) = elec(nstart:nend)+pat(kbit)*elpulse; % add the kbit-pulse
end

if flag == 1    % filter the signal
    Hf = myfilter(ptype,GSTATE.FN,roll,ord);
    elec = ifft(fft(elec).* Hf);
    delay = evaldelay(ptype,roll);
    elec=circshift(elec,-delay);
end

if strcmp(inpow,'power')
    elec = sqrt(abs(elec));
end

%------------------------------------------------
function y=pulseshape(ptype,roll,duty,Nt)

%PULSESHAPE Creates the fundamental pulse
%   Y=PULSESHAPE(PTYPE,ROLL,DUTY,NT) returns in Y a vector [2*NT,1]
%   containing the fundamental pulse whose type is defined in ptype (see
%   ELECTRICSOURCE for available types). NT is the number of points per bit,
%   DUTY the bit duty-cycle. ROLL is the pulse roll-off if used, otherwise
%   any number. The length of Y is 2*NT because the roll-off spreads the
%   pulse outside the bit time.
%   BEWARE: PULSE TYPES 'SECH' AND 'TANH' YET TO BE IMPLEMENTED

elpulse = zeros(Nt*2,1);     % elementary pulse (over two bit times
% because the roll-off spreads the pulse).
switch ptype
    case 'cosroll'

        nl = round(0.5*(1-roll)*duty*Nt);    % start index of cos roll-off
        nr = duty*Nt-nl-1;                   % end index of cos roll-off

        nmark = 1:nl;                       % indexis where the pulse is 1
        ncos  = nl:nr;                      % transition region of cos roll-off

        elpulse(Nt+nmark) = 1;
        hperiod = duty*Nt-2*nl;
        if hperiod ~= 0
            elpulse(ncos+Nt+1) = 0.5*(1+cos(pi/(hperiod)*(ncos-nl+0.5)));
        end
        elpulse(1:Nt) = flipud(elpulse(Nt+1:Nt*2)); % first half of the pulse

    case 'dirac'
        elpulse(Nt+1) = 1;
    case 'idpulse'
        nl = round(0.5*duty*Nt);
        if nl == 0
            elpulse(Nt+1) = 1;  % same as Dirac's delta. Why this? Because in this
            % way you can create a gaussian pulse, for
            % instance, by filtering this train of delta.
        else
            elpulse(Nt-nl+1:Nt+nl) = 1;
        end
    otherwise
        error('error in pulseshape.m: the pulse ptype does not exist');
end

y=elpulse;

function par=checkparam(par, format)

nch=length(format);
if isempty(par)
    % Setting default values for the parameter structure:
    par.efilt = 'bessel5';
    par.efiltbw   = 0.3;
elseif ~isfield(par,'efilt')
    par.efilt = 'bessel5';
elseif ~isfield(par,'efiltbw')
    par.efiltbw = 0.3;
end

if ~ischar(par.efilt)
    error('par.efilt must be a string');
end
if ~isnumeric(par.efiltbw)
    error('par.efilt must be a number');
end
