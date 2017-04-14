function varargout=pat_encoder(pat,modformat)
%PAT_ENCODER Symbols encoder
%   PAT=PAT_ENCODER(PAT,MODFORMAT,OPTIONS) given the pattern PAT generated 
%   by PATTERN returns the encoded pattern.
%   MODFORMAT is the modulation format and can be 'ook', 'dpsk','nf-dpsk',
%   'psbt', 'dqpsk', 'nf-dqpsk'.
%
%   Such function must be called before ELECTRICSOURCE. With 'ook'
%   it is actually unnecessary to call this function.
%
%   In case of 'dqpsk', 'nf-dqpsk', the pattern is assumed to be quaternary 
%   unless OPTIONS.binary is set to true. 
%
%   [PAT PATBIN]=PAT_DECODER(PAT,MODFORMAT,OPTIONS) returns in PATBIN also 
%   the binary representation of the M-ary pattern and thus is a matrix 
%   [NSYMB, log2(M)], being M the cardinality of the pattern.
%
%   In case of PDM coherent transmissions, PAT_ENCODER operates
%   only on the pattern associated on one polarization, since the
%   differential codings on the two polarizations are independent.
%   Thus, you have to calculate the encoded pattern for both the 
%   x and y components of the signal.
%
%   See also PATTERN, ELECTRICSOURCE, PAT_DECODER
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

if nargin > 2 && isfield(options,'binary')
    binary = options.binary;
else
    binary = false;
end

switch modformat
    case 'ook'
        % nothing to do
    case {'dpsk','dpsk-nf','psbt'} % differential coding.
        stars_t = pat2stars(pat,'dpsk');
        stars_r = -ones(size(stars_t));

        for ii=2:length(stars_r)
            stars_r(ii) = -conj(stars_t(ii)).*stars_r(ii-1);
        end

        enc = ~stars2pat(stars_r,'dpsk');

        if strcmp(modformat,'dpsk-nf')
            enc = fastshift(enc,1);
        end
    case {'dqpsk','dqpsk-nf'} % differential quadrature coding
        stars_t = pat2stars(pat,'dqpsk');
        stars_r = -ones(size(stars_t));

        for ii=2:length(stars_r)
            stars_r(ii) = -conj(stars_t(ii)).*stars_r(ii-1);
        end

        [enc encmat] = stars2pat(stars_r,'dqpsk');
        encmat=~encmat;

        if strcmp(modformat,'dqpsk-nf')
            enc = fastshift(enc,1);
            encmat = fastshift(encmat,1);
        end
    otherwise
        error('wrong modulation format in pat_encoder');
end

varargout{1}=enc;
if exist('encmat','var')
    varargout{2}=encmat;
end
