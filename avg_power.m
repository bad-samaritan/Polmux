function E=avg_power(ich,flag,fil,bw,ord)

%AVG_POWER Evaluate the average power per symbol.
%	E=AVG_POWER(ICH) returns in E the average power [mW] per symbol of 
%   channel ICH.
%
%   E=AVG_POWER(ICH,FLAG) has the optional parameter FLAG that can be 'abs' 
%   or 'norm'. FLAG='norm' allows to return the power normalized to the 
%   transmitted peak power, otherwise with FLAG='abs' (default) the power 
%   is in [mW].
%
%   E=AVG_POWER(ICH,FLAG,FIL,OBW,ORD) temporary extracts the channel using
%   the filter FIL with bandwidth OBW and order ORD (see MYFILTER).
%
%   The average power is evaluated in the frequency domain.
%
%   In the case of a unique field (see CREATE_FIELD) the function evaluates
%   the power of channel ICH (e.g. ICH may be the one centered on f1 or f2  
%   or f3 in the bottom figure) over the window D1, D2 or D3, respectively.
%
%   Spectrum
%       ^        D1                  D2                   D3
%       | <-----------------><----------------><----------------------->
%       |                    .         __      .   ______________      
%       |           /\       .       /   \     .  |              |     
%       |          /  \      .     _|     |_   .  |              |     
%       |         /    \     .    /         \  .  |              |     
%       |----------------------------------------------------------------->
%                   f1   (f1+f2)/2    f2   (f3+f2)/2    f3        frequency 
%
%
%   E contains also the contribute of GSTATE.FIELDY, if it exists.
%
%   Example: An ideal On-off keying signal with non-return to zero pulses
%       has E=0.5 with FLAG='norm'.
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
CLIGHT = CONSTANTS.CLIGHT;  % speed of light [m/s]

[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
minfreq = GSTATE.FN(2)-GSTATE.FN(1);

if (ich > GSTATE.NCH)
    error('The channel does not exist');
end
if nargin >= 2
    if strcmpi(flag,'norm')
        nrm = GSTATE.POWER(ich);
    else
        nrm = 1;
    end
else
    nrm = 1;
end
custom_filter = nargin > 2;
if nargin<5
    ord=0;
end


Nfft = length(GSTATE.FN);
if nfc ~= GSTATE.NCH            % Unique field
    maxl=max(GSTATE.LAMBDA);
    minl=min(GSTATE.LAMBDA);
    lamc = 2*maxl*minl/(maxl+minl); %central wavelength
    deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich)); % frequency spacing    
    ndfn = round(deltafn./GSTATE.SYMBOLRATE/minfreq);  % spacing in points
    nch = 1;
else                % Separate channels
    ndfn = 0; % All points are valid
    nch = ich;
end    

% ndfnl & ndfnr are used for evaluating the energy after
if (ich == 1) || (nfc == GSTATE.NCH) 
    ndfnl = Nfft/2+ndfn;
else
    deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich-1)); % left channel
    ndfnl = round(deltafn./GSTATE.SYMBOLRATE/minfreq);
    ndfnl = round((ndfn-ndfnl)*0.5);    % spacing from left chan
end
if (ich == GSTATE.NCH) || (nfc == GSTATE.NCH) 
    ndfnr = Nfft/2-ndfn;
else
    deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich+1)); % right channel
    ndfnr = round(deltafn./GSTATE.SYMBOLRATE/minfreq);        
    ndfnr = round((ndfnr-ndfn)*0.5);    % spacing from right chan
end

x = fft(GSTATE.FIELDX(:,nch));
x = fastshift(x,ndfn); % Faster than: x = x(nind);

if custom_filter
    x = x.*myfilter(fil,GSTATE.FN,bw*0.5,ord);
end


E = sum(abs(x(1:ndfnl)).^2) + sum(abs(x(Nfft-ndfnr+1:Nfft)).^2);
E = E/nrm/Nfft^2; % normalized average energy x symbol

isy = ~isempty(GSTATE.FIELDY);
if isy
    x = fft(GSTATE.FIELDY(:,nch));
    x = fastshift(x,ndfn); % Faster than: x = x(nind); 
    if custom_filter
        x = x.*myfilter(fil,GSTATE.FN,bw*0.5,ord);
    end
    Ey = sum(abs(x(1:ndfnl)).^2) + sum(abs(x(Nfft-ndfnr+1:Nfft)).^2);
    Ey = Ey/nrm/Nfft^2; % normalized average energy x symbol    
    E = E + Ey;
end
