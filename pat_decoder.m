function varargout=pat_decoder(pat,modformat,options)

%PAT_DECODER Symbols decoder
%   PAT=PAT_DECODER(PAT,MODFORMAT,OPTIONS) given the pattern PAT generated 
%   by PATTERN returns the decoded pattern.
%   MODFORMAT is the modulation format and can be 'ook', 'dpsk','nf-dpsk',
%   'psbt', 'dqpsk', 'nf-dqpsk'.
%
%   Such function must be called before functions that operates on the
%   received signal, like EVAL_EYE, BER_KL, BER_ESTIMATE, etc. With 'ook'
%   it is actually unnecessary to call this function.
%
%   In case of 'dqpsk', 'nf-dqpsk', the pattern is assumed to be quaternary 
%   unless OPTIONS.binary is set to true. 
%
%   [PAT PATBIN]=PAT_DECODER(PAT,MODFORMAT,OPTIONS) returns in PATBIN also 
%   the binary representation of the M-ary pattern and thus is a matrix 
%   [NSYMB, log2(M)], being M the cardinality of the pattern.
%
%	In case of PDM coherent transmissions, PAT_DECODER operates
%   only on the pattern associated on one polarization, since the
%   differential codings on the two polarizations are independent.
%	Thus, you have to calculate the decoded pattern for both the 
%	x and y components of the signal.
%
%   See also PATTERN, PAT_ENCODER, EVAL_EYE, BER_KL, BEST_EYE, BEST_SP,
%            BER_ESTIMATE, RECEIVER_OOK, RECEIVER_DPSK
%
%   Author: Paolo Serena, 2009
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

if nargin > 2 && isfield(options,'binary')
    binary = options.binary;
else
    binary = false;
end

switch modformat
    case 'ook'        
        % nothing to do
    case 'dpsk' % differential coding.
        stars_t = pat2stars(pat,'dpsk');
        stars_r=zeros(size(stars_t));
        stars_r=conj(stars_t).*fastshift(stars_t,1);
        pat = ~stars2pat(stars_r,'dpsk');
    case {'nf-dpsk','psbt'} % differential coding with ook receiver
        stars_t = pat2stars(pat,'nf-dpsk');
        stars_r=zeros(size(stars_t));
        stars_r=conj(stars_t).*fastshift(stars_t,1);
        pat = fastshift(~stars2pat(stars_r,'nf-dpsk'),-1);        
    case {'dqpsk','nf-dqpsk'}
        if binary
            stars_t = pat2stars(pat,'dqpsk',struct('binary',binary));
        else
            stars_t = pat2stars(pat,'dqpsk');
        end
        stars_r=zeros(size(stars_t));
        stars_r=conj(stars_t).*fastshift(stars_t,1);
        [pat patmat] = stars2pat(stars_r,'dqpsk');
        patmat=~patmat;     % both pattern have to be inverted
        pat=3-pat;
        if strcmp(modformat,'nf-dqpsk')
            pat = fastshift(pat,-1);
            patmat = fastshift(patmat,-1);
        end
    otherwise
        error('wrong modulation format in pat_decoder');
end

varargout{1}=pat;
if exist('patmat','var')
    varargout{2}=patmat;
end
