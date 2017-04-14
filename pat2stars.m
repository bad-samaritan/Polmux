function stars=pat2stars(pat,format,options)
%PAT2STARS Convert an M-ary pattern into a complex constellation
%   STARS=PAT2STARS(PAT,FORMAT) returns in STARS the values of the complex
%   constellation associated with the M-ary pattern PAT. The supported 
%   modulation formats are:
%
%       - OOK
%       - PSBT
%       - BPSK, DPSK, NF-DPSK
%       - QPSK, DQPSK, NF-DQPSK
%
%   The constellation is associated with the pattern through Gray coding,
%   thus, for example, the QPSK constellation associated with the symbols 
%   [0, 1, 2, 3] is [1, i, -i, -1].
%
%   STARS=PAT2STARS(PAT,FORMAT,OPTIONS) allows to modify the default
%   behaviour. Currently supported fields are:
%  
%   - OPTIONS.BINARY: if set to true, PAT2STARS assumes the pattern is
%   binary and thus is a matrix of size [NSYMB, log2(M)], being M the
%   cardinality of the alphabet. If M=2, this option is simply neglected.
%
%   See also STARS2PAT, PAT_DECODER
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

if nargin > 2
    if isfield(options,'binary')
        binary = options.binary;
    end
else
    binary = false;
end


stars=zeros(size(pat));
switch format
    case {'ook','psbt'}
        if max(pat)>1
            error(sprintf('pattern of %s must be binary',format));
        end
        % nothing to do
    case {'bpsk','dpsk','nf-dpsk'}
        if max(pat)>1
            error(sprintf('pattern of %s must be binary',format));
        end        
        stars(pat==0)=1;
        stars(pat==1)=-1;        
    case {'dqpsk','nf-dqpsk','qpsk'}
        if ~binary
            if max(pat)>3
                error(sprintf('pattern of %s must be quaternary',format));
            end
            stars(pat==0)=1;
            stars(pat==1)=i;
            stars(pat==3)=-1;
            stars(pat==2)=-i;
        else
            if max(max(pat))>1 || size(pat,2) ~= 2
                error('pattern must be a binary matrix with size [Nsymb,2]');
            end
            stars=zeros(length(pat),1);
            
            for ii=1:length(stars)
                if pat(ii,:)==[0 0] stars(ii)=1;end
                if pat(ii,:)==[0 1] stars(ii)=i;end
                if pat(ii,:)==[1 1] stars(ii)=-1;end
                if pat(ii,:)==[1 0] stars(ii)=-i;end                
            end
        end        
    otherwise
        error('unknown modulation format')        
end
