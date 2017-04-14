function Eout=linear_modulator( Ein, modsig, exr)
%LINEAR_MODULATOR modulate the optical field with a linear modulator
%   E=LINEAR_MODULATOR(E,MODSIG) modulates the optical field E using a 
%   linear modulator, i.e. the modulating electrical signal is impressed 
%   exactly on the optical signal with an infinite extinction ratio.
%   The parameter MODSIG is the electrical driving signal produced by 
%   ELECTRICSOURCE. 
%   E=LINEAR_MODULATOR(E,MODSIG,EXRATIO) adds a finite extinction ratio 
%   EXRATIO [dB] to the signal.
%   
%   See also LASERSOURCE, ELECTRICSOURCE, MZ_MODULATOR
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


exratio = inf;

if exist('exratio','var');
    exratio = exr;
end

exr_lin = 10^(-exratio/10);

Eout = Ein.*((1-sqrt(exr_lin)).*modsig + sqrt(exr_lin));
