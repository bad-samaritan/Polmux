function varargout=stars2pat(stars,format)
%STARS2PAT Convert a complex constellation into a pattern 
%   PAT=STARS2PAT(STARS,FORMAT) returns in PAT the M-ary pattern associated 
%   with the complex constellation. The supported modulation formats are:
%
%       - OOK
%       - PSBT
%       - BPSK, DPSK, NF-DPSK
%       - QPSK, DQPSK, NF-DQPSK
%
%   The association rule between constellation and pattern is the same as
%   in PAT2STARS.
%
%   [PAT BPAT]=STARS2PAT(STARS,FORMAT) returns also the binary equivalent
%   of the M-ary pattern.
%
%   See also PAT2STARS, PAT_DECODER
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
pat = zeros(size(stars));

switch format
    case {'ook','psbt'}
        % nothing to do
    case {'bpsk','dpsk','nf-dpsk'} 
        pat(stars==1)=0;
        pat(stars==-1)=1;        
    case {'dqpsk','nf-dqpsk','qpsk'}          
        pat(stars==1)=0;
        pat(stars==i)=1;
        pat(stars==-1)=3;
        pat(stars==-i)=2;
        
        patmat=zeros(length(pat),2);
        
        patmat(stars== 1,1)=0; patmat(stars== 1,2)=0;
        patmat(stars== i,1)=0; patmat(stars== i,2)=1;
        patmat(stars==-1,1)=1; patmat(stars==-1,2)=1;
        patmat(stars==-i,1)=1; patmat(stars==-i,2)=0;        
    
        varargout{2}=patmat;
    otherwise
        error('unknown modulation format')               
end

varargout{1}=pat;
