function Eout=lpfilter(Ein,ftype,bw,ord)

%LPFILTER Filtering with a lowpass filter.
%   EOUT=LPFILTER(EIN,FTYPE,BW) filters the input field EIN with the filter 
%   FTYPE (see MYFILTER) having 3-dB bandwidth BW (normalized to the 
%   symbolrate). The resulting filtered field is associated to EOUT 
%
%   EOUT=LPFILTER(EIN,FTYPE,BW,ORD) use the additional parameter ORD for special
%   filter. E.g. ORD is the supergauss order for the supergaussian filter 
%   (see MYFILTER)
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

global GSTATE

if nargin < 4
    ord = 0;
end


Hf = myfilter(ftype,GSTATE.FN,bw,ord); % Remember that the in the lowpass
%   equivalent domain, the 3 dB bandwidth goes from -bw/2 to + bw/2
delay = evaldelay(ftype,bw);
shift=round(delay*GSTATE.NT);

Eout = ifft(fft(Ein) .* Hf);
Eout = circshift(Eout,-shift);
