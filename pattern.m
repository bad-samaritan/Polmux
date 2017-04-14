function varargout=pattern(ptype,nseed,options)

%PATTERN Create the sequence pattern with rules.
%   PAT=PATTERN(PTYPE,NSEED,OPTIONS) returns in PAT a sequence of Nsymb 
%   integers representing the symbol-pattern for a digital modulation. 
%   Nsymb=GSTATE.NSYMB see RESET_ALL. In absence of OPTIONS the pattern is 
%   a bit-pattern, i.e. a vector of 0 and 1. 
%
%   PTYPE is the type of the pattern and can be one of the following:
%
%       'debruijn': creates a De Bruijn sequence (DBS) (each subsequence of
%               length log2(Nsymb) appears exactly once in a DBS) [1-2]. In 
%				the binary case, a DBS is a pseudo random binary sequence 
%			    (PRBS) with an additional zero added to to the longest 
%				sequence of 0. 
%
%               NSEED is the DBS seed. 0 <= NSEED < Nsymb/4 yields a unique 
%				DBS, i.e. it is not possible to obtain the same DBS with a 
%				different NSEED, neither with a circular shift. 
%               For NSEED >= Nsymb/4 the DBS is not unique, hence PATTERN
%               generates a circular random delayed version of a DBS 
%               with NSEED < Nsymb/4. 
%   	        NSEED must be < Nsymb/4*(Nsymb-1).
%               It is not possible to have the same sequence for different 
%               NSEED.
%
%       'random': creates a uniform distributed random sequence.
%
%       <sequence of numbers>: E.g. '20104101' creates a periodic 
%                 repetition of the sequence up to length Nsymb, and 
%                 truncates when necessary.  
%                 For instance, with Nsymb=16 returns 2010410120104101. 
%                 The sequence can be a string or a vector of double, e.g.
%                 '20104101' or [2 0 1 0 4 1 0 1]
%
%       <file>: reads the pattern from 'file' using LOAD. 'file' contains 
%               the pattern for the channel.
%               E.g. Nsymb=8, with patterns '30101210', 'file' can have any 
%               form, e.g.:
%
%                   3 0 1 0    or        3 0 1 0 1 2 1 0 
%                   1 2 1 0              
%
%
%   OPTIONS is an optional parameter containing:
%
%       OPTIONS.alphabet: is the alphabet of the pattern. 
%           E.g. PAT=PATTERN('debruijn',0,OPTIONS) with OPTIONS.alphabet=4 
%           returns:
%
%           PAT=[1  1  2  3  0  3  1  3  3  2  2  1  0  2  0  0]
%
%           i.e. a DBS sequence with symbols (0,1,2,3) containing all 
%           couples of symbols exactly once.
%
%           E.g. PAT=PATTERN('random',OPTIONS) with OPTIONS.alphabet=8 
%           may return:
%
%           PAT=[6  3  7  6  6  2  3  7  3  2  5  4  4  2  0  2]
%
%
%   [PAT,BMAT] = PATTERN(PTYPE,NSEED,OPTIONS) returns in BMAT the binary
%   representation of PAT. BMAT is a matrix of size 
%	[Nsymb,ceil(log2(OPTIONS.alphabet))].
%
%   Example: For a QPSK modulation, the two columns of BMAT represent the 
%   in-phase and quadrature component of a pseudo random quaternary 
%   sequence (PRQS), [2].
%
%   Note: PAT and BMAT are of type double.
%
%   See also ELECTRICSOURCE, PAT_DECODER
%
%
%   [1] F. S. Annexstein, "Generating De Bruijn Sequences: An Efficient
%   Implementation," IEEE Transaction on Computers, vol. 48, no. 2, 
%   pp. 198-200, Feb. 1997.
%
%   [2] D. v. d. Borne et al, "Bit pattern dependence in optical DQPSK
%   modulation," Electronic Letters, vol. 43, no. 22, Oct. 2007
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


global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m


fromfile = 0;
     
if strcmp(ptype,'debruijn') % PSEUDO-RANDOM BINARY SEQUENCE
    
    if nargin == 1
        error('Missing nseed');
    else
        if nargin == 3
            nfnames = fieldnames(options);
            if any(strcmp(nfnames,'alphabet'))
                q = log2(options.alphabet);
                if rem(q,1)
                    error(['De Bruijn sequence can be implemented only for ',...
                        ' power of 2 alphabets']);
                end
            else
                q = 1;
            end
        else
            q = 1; % default: binary
        end
        maxseed = GSTATE.NSYMB*(GSTATE.NSYMB-1)/4; % largest seed
        if nseed(1) > maxseed, error('nseed must be <= %d',maxseed); end
        lg2 = log2(GSTATE.NSYMB);
        if rem(lg2,q) 
            error(['The De Bruijn sequence does not exist! log2(number ',... 
                'of symbols)=%d is not a multiple of log2(alphabet length)=%d.'],...
                lg2,q);
        end
        [pat,tmat] = debruijn_seq(log2(GSTATE.NSYMB)/q,nseed(1),q);
        varargout{2} = tmat';
    end
    
    
elseif strcmp(ptype,'random') % RANDOM PATTERN  
    if nargin >= 2
        if nargin == 2, options = nseed;end
        nfnames = fieldnames(options);
        if any(strcmp(nfnames,'alphabet'))
            qq = options.alphabet;
        else
            qq = 2;
        end
    else
        qq = 2;
    end
    pat = floor(qq.*rand(1,GSTATE.NSYMB)); 
    varargout{2} = mydec2bin(pat,log2(qq));
    
elseif isnumeric(ptype)     % VECTOR OF SYMBOLS
    if nargin >= 2
        if nargin == 2, options = nseed;end
        nfnames = fieldnames(options);
        if any(strcmp(nfnames,'alphabet'))
            q = log2(options.alphabet);
            if rem(q,1)
                error(['Sequence can be implemented only for ',...
                    ' power of 2 alphabets']);
            end
        else
            q = 1;
        end
    else
        q = 1; % default: binary
    end
    pat = myseq(ptype,GSTATE.NSYMB); % periodically repeat the pattern
    varargout{2} = mydec2bin(pat,q);
    
elseif (min(double(ptype-48)) >= 0 && max(double(ptype-48)) <= 9)   % STRING
    if nargin >= 2
        if nargin == 2, options = nseed;end
        nfnames = fieldnames(options);
        if any(strcmp(nfnames,'alphabet'))
            q = log2(options.alphabet);
            if rem(q,1)
                error(['Sequence can be implemented only for ',...
                    ' power of 2 alphabets']);
            end
        else
            q = 1;
        end
    else
        q = 1; % default: binary
    end   
    ptype2 = double(ptype-48);
    pat = myseq(ptype2,GSTATE.NSYMB); % periodically repeat the pattern
    varargout{2} = mydec2bin(pat,q);

elseif exist(ptype,'file') == 2
    
    pat = loadpattern(ptype,GSTATE.NSYMB);
    fromfile = 1;    
    varargout{2} = mydec2bin(pat);
    
else
    
    error('String ptype unknown');
end
varargout{1} = pat(:);

if GSTATE.PRINT

    %%%%% PRINT

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===             pattern              ===\n');
    fprintf(fid,'========================================\n\n');
    if strcmp(ptype,'debruijn')
        fprintf(fid,'%s pattern with seed %d\n',ptype,nseed(1));
    elseif strcmp(ptype,'random')
        fprintf(fid,'%s pattern\n',ptype);
    elseif fromfile == 1 
        fprintf(fid,'user-defined pattern from file %s\n',ptype);
    else
        fprintf(fid,'Periodic pattern:\n');
        if isnumeric(ptype), ptype=char(ptype+48);end
        fprintf(fid,'...%s...',ptype);
    end
    if exist('options','var') && isfield(options,'alphabet')
        fprintf(fid,'Alphabet: %d\n',options.alphabet);
    end
    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT


%--------------------------------------------------------------------------
function y=myseq(locpat,Nsymb)

% creates a periodic repetition of the sequence locpat up to Nsymb, and
% truncates if necessary.
%
% E.g.  ptype = [0 1 0] and Nsymb = 8 -> y=[0 1 0 0 1 0 0 1]
%       ptype = [0 1 0 0 1 0] and Nsymb = 4 -> y=[0 1 0 0]
%
%   Author: Paolo Serena, 2009
%   University of Parma, Italy

Lp = length(locpat);
nrem = rem(Nsymb,Lp);
if Lp > Nsymb
    y = locpat(1:Nsymb);
elseif nrem ~= 0    
    temp1 = repmat(locpat,1,floor(Nsymb/Lp));   % replicate    
    temp2 = locpat(1:nrem);
    y = [temp1,temp2];
else
    y = repmat(locpat,1,Nsymb/Lp);
end


%--------------------------------------------------------------------------
function y=loadpattern(ptype,Nsymb)

% load the pattern from file ptype using the rules of function load
%
%   Author: Paolo Serena, 2007
%   University of Parma, Italy


locpat = load(ptype);
[nr,nc] = size(locpat);
fprintf('%d %d %d\n',nr,nc,Nsymb)
if (nr*nc) ~= (Nsymb)
    error(['File ptype ',...
        'contains a wrong number of symbols']);
end
y = reshape(locpat',1,Nsymb);

%--------------------------------------------------------------------------
function [y,tmat]=debruijn_seq(n,seed,q)

%DEBRUIJN_SEQ De Bruijn sequence generator
%   Y=DEBRUIJN_SEQ(N,SEED) generates a De Bruijn sequence of 2^N bits.
%   A De Bruijn sequence contains all patterns of N bits exactly once.
%
%   SEED is the random seed (integer). For SEED > 2^(N-2) the resulting 
%   sequence is a pseudo-random shifted copy of a sequence with 
%   SEED <= 2^(N-2).
%
%   Reference:
%
%   [1] F. S. Annexstein, "Generating De Bruijn Sequences: An Efficient
%   Implementation," IEEE Transaction on Computers, vol. 48, no. 2, Feb.
%   1997.

ns = q*n; % log2 length of the DeBruijn sequence
N = 2^(ns-2);
nseed = mod(seed,N);
x = double(dec2bin(nseed,ns-2)-48); % convert to pattern
y = double(generate_debruijn(ns,x));

if q > 1 % multilevel DeBruijn
    mvm = floor(ns/q*((1:q)-1));
    tmat = zeros(q,2^ns); % columns: binary DeBruijn
    for kk=1:q
        tmat(kk,:) = fastshift(y,mvm(kk));
    end
    y = 2.^(q-1:-1:0) * tmat ; % bin2dec conversion
else
    tmat = [];
end
if seed > N % no more seeds -> apply a random delay shift 
    nseed2 = ceil((seed-N+1)/N);
    nshift = mod(97*nseed2,2^ns-1)+1; % congruent random generator
    y = fastshift(y,-nshift); % shift on the right
    tmat = fastshift(tmat.',-nshift).';
end

%--------------------------------------------------------------------------
function y=generate_debruijn(n,x)

if n == 2
    y=[1,1,0,0] == 1;
elseif n == 3
    if x(1) == 0
        y = [1,0,1,1,1,0,0,0] == 1;
    else
        y = [1,1,1,0,1,0,0,0] == 1;
    end
else
    y=next_debruijn(generate_debruijn(n-1,x),2^(n-2)+(-1)^x(n-3),x(n-2));
end

%--------------------------------------------------------------------------
function y=next_debruijn(w,i,k)

C = xor_iscan(w);
Cbar = not(C);
part1 = C(1:i-k);
part2 = Cbar(i+k:end);
part3 = Cbar(1:i-1+k);
part4 = C(i-k+1:end);
y = [part1,part2,part3,part4];

%--------------------------------------------------------------------------
function y=xor_iscan(x)

% Being x=[x1,x2,...,xn] and y=[y1,y2,...,yn] it is:
%       yk = x1 xor x2 xor ... xk

z = cumsum(x);
y = rem(z,2) == 1;

%--------------------------------------------------------------------------
function y=mydec2bin(d,q)

[f,e]=log2(max(d)); %#ok How many digits do we need to represent the numbers?
y=rem(floor(d(:)*pow2(1-max(1,e):0)),2);
cols = size(y,2);
if q > cols
    y(:,q-cols+(1:cols)) = y(:,1:cols);
    y(:,1:(q-cols)) = zeros(size(y,1),q-cols);
end
