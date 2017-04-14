function E=qi_modulator(E, modsig_i, modsig_q, options)
%QI_MODULATOR Modulate the optical field using a QI Mach-Zehnder modulator
%   E=QI_MODULATOR(E,MODSIG_I,MODSIG_Q, OPTIONS) modulates the optical field
%   using a QI Mach-Zehnder interferometer [1].
%   The parameters MODSIG_I and MODSIG_Q are the in-phase and in-quadrature
%   electrical driving signals produced by ELECTRICSOURCE. OPTIONS is an 
%   optional structure whose fields can be:
%
%       - iqratio : change the power ratio between I and Q arms 
%         (default = 0)
%       - biasc : change the bias between I and Q arms (default = 0)
%       - exratio : extinction ratio of the two nested modulators [dB] 
%         (default = [inf inf])
%       - bias  : bias of the two nested modulators (default = [0 0])
%       - amplitude: Vpi of the two nested modulators (default = [1 1])
%       - nochirp: reduce effect of chirp due to finite exratio of 
%         two nested modulators [2] (default = [false false])         
%   
%
%   See also PATTERN, ELECTRICSOURCE, LASERSOURCE, MZ_MODULATOR
%
%   [1] R. A. Griffin, "Integrated DQPSK transmitters," in Proc. OFC 2005,
%   2005 Anaheim, CA, USA
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

iqratio = 0;
biasc     = 0;
bias      = [0 0];
amplitude = [1 1];
exratio     = [inf inf];
nochirp   = [false false];

if exist('options','var')
    checkfields(options,{'iqratio','biasc','bias','amplitude','exratio','nochirp'});
    if isfield(options,'iqratio');
        iqratio = options.iqratio;
    end    
    if isfield(options,'biasc');
        biasc = options.biasc;
    end
    if isfield(options,'bias');
        bias = options.bias;
        if length(bias) ~= 2
            error('bias must be a vector of length2')
        end
    end    
    if isfield(options,'amplitude');
        amplitude = options.amplitude;
        if length(amplitude) ~= 2
            error('amplitude must be a vector of length2')
        end        
    end
    if isfield(options,'exratio');
        exratio = options.exratio;
        if length(exratio) ~= 2
            error('exratio must be a vector of length2')
        end        
    end  
    if isfield(options,'nochirp');
        nochirp = options.nochirp;
        if length(nochirp) ~= 2
            error('nochirp must be a vector of length2')
        end        
    end      
end

iqratio = 10^(iqratio/20);
sr = iqratio/(1+iqratio);

%   The main function:
Ei = E * sr;         % In Phase Field
Eq = E * (1-sr);     % In Quadrature Field

Ei = mz_modulator(Ei,modsig_i,struct('exratio',exratio(1),'bias',...
    bias(1),'amplitude',amplitude(1),'nochirp',nochirp(1)));
Eq = mz_modulator(Eq,modsig_q,struct('exratio',exratio(2),'bias',...
    bias(2),'amplitude',amplitude(2),'nochirp',nochirp(2)));

E = sqrt(2)*(Ei+Eq*exp(i*(pi/2 + biasc)));  
