function [cond,avgber,nruns,stdber]=ber_estimate(pat_hat,pat,x,nind)

%BER_ESTIMATE Bit-error rate estimate by Monte Carlo simulation
%   [COND,AVGBER,NRUNS,STDBER]=BER_ESTIMATE(PAT_HAT,PAT,X) estimates the 
%   average bit error rate (AVGBER) by standard Monte Carlo (MC) simulation. 
%
%   NRUNS is the number of samples used in the estimation. STDBER is the
%   standard deviation of the measured BER. COND is true during the MC  
%   cycle, and false when the convergence stop criterion has been reached.
%   PAT is the bit pattern, decoded with pat_decoder.m if necessary, while
%   PAT_HAT is the received estimated pattern given by SAMP2PAT.
%
%   X is a struct with one or both of the following fields:
%
%   X.stop = vector [t1 t2]. t1 is the relative BER accuracy, while t2 the  
%       Gaussian confidence of the accuracy. In such a case the MC 
%       simulation tests a sufficient number of runs as soon as the
%       accuracy is reached. For instance, with t1=0.01 and t2=95, the 
%       simulation ends as soon as the ratio  error/AVGBER < 0.01, being 
%       error the estimated error on AVGBER with a Gaussian confidence of 
%       95%. Beware that the concept of confidence works for Gaussian 
%       distributed random variables, while in the general case with few
%	samples it is just a rule of thumb.  
%
%       A Gaussian confidence of 68% means that a Gaussian distributed BER 
%       is within +/- STDBER with probability 68%.
%
%   X.nmin = minimum number of errors counted. In absence of X.stop this is
%       also the overall number of errors counted. If not specified, the 
%       default value is 1.
%
%   [COND,AVGBER,NRUNS,STDBER]=BER_ESTIMATE(PAT_HAT,PAT,X,NIND) runs a 
%   vectorial MC. In this case X must have the additional field:
%
%   X.dim = dimension of AVGBER (external cycle dimension).
%
%   COND,AVGBER,NRUNS,STDBER are vectors of size X.dim, and they are 
%   currently evaluated in position (or index) NIND.
%
%   Such options is useful for MC estimation inside a loop and allows to 
%   evaluate all the entries of AVGBER with the same accuracy given by 
%   X.stop and X.nmin. COND is a vector of logical. 
%
%   Example
%
%   MC estimation for a system by testing a vector of post-compensating 
%   fibers placed at the end of the link. See the examples for more 
%   information.
%
%   See also EVAL_EYE, BER_KL, PDF2BER, PDF_ESTIMATE, PAT_DECODER
%
%   Authors: Paolo Serena, 2009
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

%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION

fnames=fieldnames(x);

if ~any(strcmp(fnames,'nmin')), x.nmin = 1; end
if ~any(strcmp(fnames,'dim'))
    if nargin < 4
        x.dim = 1;
        nind = 1;
    else
        error('missing variable nind in ber_estimate');
    end
end
if any(strcmp(fnames,'stop'))
    if (x.stop(2) > 100) || (x.stop(2) < 0)
        error('the Gaussian confidence must be < 100 and > 0');
    end
    isepsilon=true;
else
    isepsilon=false;
end

[cond,avgber,nruns,stdber]=mc_run(pat_hat,pat,x,isepsilon,nind);    


%--------------------------------------------------------------------------
function [cond2,avgber2,n2,stdber]=mc_run(pat_hat,pat,x,isepsilon,nind)

%%%%%%%%%%%%%%%%%%%%%% MC SIMULATION
% The simulation collects M samples each iteration, hence it is a block
% iteration.

persistent n avgber varber first epsilon cond; % keep memory

M = size(pat,1)*size(pat,2); % block dimension
if isempty(first) 
    n = ones(1,x.dim);  % first time only: initialization
    avgber = zeros(1,x.dim);
    varber = zeros(1,x.dim);
    first = 1; % remember for next time (persistent variable)
    cond = true(1,x.dim);  % vector of trues
    if isepsilon
            epsilon = sqrt(2)*erfcinv(1-x.stop(2)/100);
    end
end
nnew = n(nind)*M;

err = sum(sum(pat ~= pat_hat));

N = (n-1)*M; % cumulated number of samples
varerr = (err - err^2/M)/(M-1); % variance of the current block of M samples
avgerr = err/M; % mean of the current block
% varber(nind) is the cumulative variance
varber(nind) = ((N(nind)-1)*varber(nind) + (M-1)*varerr + ...
    N(nind)*M/(N(nind)+M)*(avgber(nind)-avgerr)^2)/(N(nind)+M-1);  
avgber(nind) = ((n(nind)-1)*avgber(nind) + avgerr)/n(nind); % cumulative mean
stdber = sqrt(varber./(N+M)); % standard deviation. 
if isepsilon       
    % Note: the N+M is the MC decreasing variance factor
    if (epsilon*stdber(nind) < x.stop(1)*avgber(nind)) && (avgber(nind)*nnew >= x.nmin)
        cond(nind) = false; % stop MC
        if ~any(cond), first = []; end
    end
else
    if avgber(nind)*nnew > x.nmin 
        cond(nind) = false;
        if ~any(cond), first = []; end
    end
end
n(nind) = n(nind) + 1;   
n2 = (n-1)*M;
avgber2 = avgber;
cond2 = cond;
