function pat_rx = samp2pat(x,s,outvalue)

%SAMP2PAT Convert received samples into a pattern
%   PAT_RX = SAMP2PAT(X,S,OUTVALUE) takes the samples OUTVALUE returned by 
%   EVAL_EYE or DSP4COHDEC and converts them into a receiver pattern
%   PAT_RX. Currently supported modulation formats are 'ook', 'dpsk',
%   'dqpsk', 'qpsk'.
%   X is a structure that describes the receiver, while S is a structure
%   used to specify how to convert the samples into a pattern. Its fields 
%   are: 
%       S.alphabet = vector of the alphabet of the desired output pattern 
%       S.thr      = vector of threshold values to apply in order to decide 
%                    how to map a sample into an alphabet symbol. S.thr
%                    length must be equal to S.alphabet length minus one
%   
%   See also: EVAL_EYE, DSP4COHDEC, BER_ESTIMATE
%

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
% Initial checks

fnames = fieldnames(s);

if any(strcmp(fnames,'thr')), isthr=false;end; 
s.logic=sort(s.logic,'ascend'); % sort symbols    
s.La = length(s.logic);
if ~isthr 
    if length(s.thr) ~= s.La-1
        error('The number of thresholds must be the alphabet length - 1');
    else
        thr = sort(s.thr,'ascend');
    end
end

% Convert samples to pattern

switch x.rec
    case {'ook','dpsk','dqpsk'}
        pat_rx = zeros(size(outvalue));
        pat_rx(outvalue >  s.thr(1))=s.logic(2);
        pat_rx(outvalue <= s.thr(1))=s.logic(1);
    case 'coherent'
        second_bit = outvalue > 0;
        first_bit = abs(outvalue) <= pi/2;
        if size(outvalue,2) == 1
            pat_rx = [first_bit second_bit];
        else
            pat_rx = [first_bit(:,1) second_bit(:,1)  first_bit(:,2)  second_bit(:,2) ];
        end
    otherwise
        error('Wrong modulation format\n');
end        
