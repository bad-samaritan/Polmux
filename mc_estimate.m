function [cond,out]=mc_estimate(s,x,nind)

%MC_ESTIMATE Monte Carlo estimation of a random variable mean and variance
%   [COND,OUT]=MC_ESTIMATE(S,X) estimates the average value and variance of 
%   a random variable by standard Monte Carlo (MC) simulation. 
%
%   S is the vector of random samples. 
%   X is a struct with one or both of the following fields:
%
%   X.stop = vector [t1 t2]. t1 is the relative accuracy of the estimator, 
%       while t2 the Gaussian confidence of the accuracy. In such a case the 
%       MC simulation tests a sufficient number of runs as soon as the
%       accuracy is reached. For instance, with t1=0.01 and t2=95, the 
%       simulation ends as soon as the ratio  OUT.stdmean/OUT.mean < 0.01, 
%       if X.method = 'mean'. Otherwise the simulation ends with 
%       (O2-O1)/OUT.var < 0.01 if X.method = 'var' (see OUT.varlim next).
%       In the previous example the Gaussian confidence was 95%. Beware that 
%       the concept of confidence works for Gaussian distributed random 
%       variables, while in the general case with few samples it is just a 
%	rule of thumb. 
%
%       A Gaussian confidence of 68% means that the exact mean of a Gaussian 
%	distributed S is within +/- OUT.stdmean of OUT.mean with probability 
%	68%.
%
%   X.nmin = minimum number of samples for applying X.stop. 
%       In absence of X.stop this is also the overall number of samples
%       tested.
%    
%   X.method = 'mean' means that X.stop is applied on the estimated average 
%	value. 'var' means that X.stop is applied to the variance. 
%	Default is 'mean'.
%
%   On output the function returns COND which is true during the MC  
%   cycle, and false when the convergence stop criterion has been reached.
%
%   OUT is a struct containing the results, i.e.:
%
%   OUT.mean = estimated average value of X;
%   OUT.var = estimated variance of X;
%   OUT.nruns = overall number of samples tested;
%   OUT.stdmean = standard deviation of OUT.mean. Note: this is not
%       sqrt(OUT.var), but the standard deviation of the average value.
%       Hence, for increasing OUT.nruns, OUT.stdmean decreases making
%       OUT.mean a good measure of the exact average value;
%   OUT.varlim = [O1 O2]. The unknown exact variance satisfies:
%           
%           O1 < exact variance < O2
%
%       with confidence X.stop(2). Note: O1 and O2 are functions of 
%       chi square percentiles. Since for OUT.nruns > 50 they can be
%       safely approximated with Gaussian percentiles we adopted such
%       approximation (see [1]).
%
%
%   [COND,OUT]=MC_ESTIMATE(S,X,NIND) runs a vectorial MC. In this case X 
%   must have the additional field:
%
%   X.dim = dimension of OUT.mean (external cycle dimension).
%
%   COND,OUT.mean,OUT.nruns,OUT.var,OUT.stdmean are vectors of size X.dim, 
%   and they are currently evaluated in position (or index) NIND. S can be
%   a matrix or a vector. In the matrix case, all samples of column NIND 
%   are used to update the estimation of OUT at index NIND, while in the
%   vector case the entire vector (row or column) is used to update index 
%   NIND of OUT.
%
%   Such options is useful for MC estimation inside a loop and allows to 
%   evaluate all the entries of OUT.mean with the same accuracy given by 
%   X.stop and X.nmin. COND is a vector of logical. 
%
%   MC_ESTIMATE is a generalization of BER_ESTIMATE.
%
%   See also BER_ESTIMATE, EVAL_EYE, BER_KL, PDF2BER, PAT_DECODER
%
%   [1] A. Papoulis, Probability, Random Variables, and Stochastic
%   Processes", McGraw Hill. 
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
DEFMIN = 50;    % default minimum number of samples

fnames=fieldnames(x);

if ~any(strcmp(fnames,'nmin')), x.nmin = DEFMIN; end
if any(strcmp(fnames,'method'))
    if ~strcmp(x.method,'mean') && ~strcmp(x.method,'var') 
        error('field method must be ''mean'' or ''var''');
    end
else
    x.method = 'mean'; 
end

if ~any(strcmp(fnames,'dim'))
    if nargin < 3
        x.dim = 1;
        nind = 1;
    else
        error('missing variable nind mc_estimate');
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

[cond,out]=complete_mc(s,x,isepsilon,nind);    


%--------------------------------------------------------------------------
function [cond2,out]=complete_mc(s,x,isepsilon,nind)

%%%%%%%%%%%%%%%%%%%%%% MC SIMULATION
% The simulation collects M samples each iteration, hence it is a block
% iteration.
%
% Ref: [1]: A. Papoulis, Probability, Random Variables, and Stochastic
%   Processes", McGraw Hill, 3rd edition.

persistent n avgout varout first epsilon cond varlim; % keep memory

if min(size(s)) == 1
    s = s(:);
    nind2 = 1;
else
    nind2=nind;
end
M = size(s,1); % block dimension
if isempty(first) 
    n = ones(1,x.dim);  % first time only: initialization
    avgout = zeros(1,x.dim);
    varout = zeros(1,x.dim);
    varlim = zeros(2,x.dim);
    first = 1; % remember for next time (persistent variable)
    cond = true(1,x.dim);  % vector of trues
    if isepsilon
        epsilon(1) = sqrt(2)*erfcinv(1-x.stop(2)/100);
        epsilon(2) = sqrt(2)*erfcinv(1+x.stop(2)/100);
    else
        epsilon = [0 0]; % not used in this case
    end
end
runs = n*M;

N = (n-1)*M; % cumulated number of samples
if (M == 1) && (n(nind) == 1)
    varout(nind) = 0; % just one sample
    avgout(nind) = s;
else
    varblk = var(s); % variance of the current block of M samples
    avgblk = mean(s); % mean of the current block
    % varout(nind) is the cumulative variance
    varout(nind) = ((N(nind)-1)*varout(nind) + (M-1)*varblk(nind2) + ...
        N(nind)*M/(N(nind)+M)*(avgout(nind)-avgblk(nind2))^2)/(N(nind)+M-1);
    avgout(nind) = ((n(nind)-1)*avgout(nind) + avgblk(nind2))/n(nind); % cumulative mean
end
out.stdmean = sqrt(varout./(N+M)); % standard deviation. 
x21mdh = 0.5*(epsilon(1) + sqrt(2*(N(nind)+M)-3))^2; % [1], p. 253 below the table
x2dh = 0.5*(epsilon(2) + sqrt(2*(N(nind)+M)-3))^2;
varlim(1,nind) = (N(nind)+M-1)*varout(nind)/x21mdh; % see [1], p. 253.
varlim(2,nind) = (N(nind)+M-1)*varout(nind)/x2dh;

if isepsilon       
    % Note: the N+M is the MC decreasing variance factor
    if strcmp(x.method,'mean') % average test
	    absavg = abs(avgout(nind));
        if (epsilon(1)*out.stdmean(nind) < x.stop(1)*absavg) && ...
                (runs(nind) >= x.nmin)
            cond(nind) = false; % stop MC
            if ~any(cond), first = []; end
        end
    else                       % variance test
        if (varlim(2,nind)-varlim(1,nind))/varout(nind) < x.stop(1) &&...
                (runs(nind) >= x.nmin)
            cond(nind) = false; % stop MC
            if ~any(cond), first = []; end
        end
    end        
else
    if runs(nind) > x.nmin 
        cond(nind) = false;
        if ~any(cond), first = []; end
    end
end
n(nind) = n(nind) + 1;   
out.nruns = (n-1)*M;
out.mean = avgout;
out.var = varout;
out.varlim = varlim;
cond2 = cond;
