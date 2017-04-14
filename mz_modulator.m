function Eout=mz_modulator( Ein, modsig, options)
%MZ_MODULATOR modulate the optical field with a Mach-Zehnder Interferometer
%   E=MZ_MODULATOR(E,MODSIG, OPTIONS) modulates the optical field E using a 
%   Mach-Zehnder interferometer [1]. 
%   The parameter MODSIG is the electrical driving signal produced by 
%   ELECTRICSOURCE. OPTIONS is an optional structure whose fields can be:
%
%       - exratio : extinction ratio [dB] (default = inf)
%       - bias  : bias of the modulator (default = 0)
%       - amplitude: Vpi of the modulator (default = 1)
%       - nochirp: reduce effect of chirp due to finite exratio [2]
%         (default = false)
%   
%   See also LASERSOURCE, ELECTRICSOURCE, LINEAR_MODULATOR, PHASE_MODULATOR
%   
%   [1] S. Walklin and J. Conradi, "Effect of Mach-Zehnder modulator DC 
%   extinction ratio on residual chirp-induced dispersion in 10-Gb/s binary
%   and AM-PSK duobinary lightwave systems," IEEE Photon. Technol. Lett., 
%   vol. 9, pp. 1400-1402, Oct. 1997.
%   [2] K. Hoon and A. H. Gnauck, "Chirp characteristics of dual-drive 
%   Mach-Zehnder modulator with a finite DC extinction ratio," IEEE Photon. 
%   Technol. Lett., vol. 14, pp. 298-300, Mar. 2002.
%
%   Author: Massimiliano Salsi, 2009
%   Contributed by: Marco Bertolini, 2009
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


bias      = 0;
amplitude = 1;
exratio     = inf;
nochirp   = false;

if exist('options','var')
    checkfields(options,{'bias','amplitude','exratio','nochirp'});
    if isfield(options,'bias');
        bias = options.bias;
    end
    if isfield(options,'amplitude');
        amplitude = options.amplitude;
    end
    if isfield(options,'exratio');
        exratio = options.exratio;
    end  
    if isfield(options,'nochirp');
        nochirp = options.nochirp;
    end      
end

exr_lin = 10^(-exratio/10);
gamma = (1-sqrt(exr_lin))/(sqrt(exr_lin)+1);

if nochirp
    Phi_U = (modsig *          pi/(1+gamma^2) * amplitude + bias);
    Phi_L = (modsig * -gamma^2*pi/(1+gamma^2) * amplitude + bias);
else
    Phi_U = (modsig *  pi/2 * amplitude + bias);
    Phi_L = (modsig * -pi/2 * amplitude + bias);
end

% Any phase shift due to the interferometric structure is neglected for the
% sake of simplicity
Eout  = j*Ein.*(fastexp(Phi_L)-gamma*fastexp(Phi_U))/(1+gamma);
