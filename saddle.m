function [u0,nsw] = saddle(sgn,xi,dx,dsg2,qsg2,qsg4,vars,ld,ldex,...
    ld2,b2,pol2)

%SADDLE Evaluate the MGF saddle point
%   [U0,NSW] = SADDLE(SGN,XI,DX,DSG2,QSG2,QSG4,VARS,LD,LDEX,LD2,B2,POL2,NSYMB)
%   returns the real saddlepoint U0 for all bits NSYMB of the moment 
%   generating function (MGF) of the sampled signal. NSW is a control 
%   parameter equal to:
%
%   NSW=0 if the saddlepoint has been evaluated with the desired accuracy
%   NSW=1 if the function did not found the saddlepoint
%   NSW=2 if the function found a saddlepoint with low accuracy
%
%   SGN is a flag equal to 1 for polar bit-formats, while for unipolar 
%   formats it takes 1 for spaces and -1 for marks. 
%   XI is (22) of [1], DX is XI-eta, being eta eq.(19) of [1]. 
%   DSG2,QSG2,QSG4 are respectively 2*sigma^2,4*sigma^2,4*sigma^4, being 
%   sigma^2 the noise variance of the sampled signal.
%
%   VARS is (20) of [1]. LD and LD2 are the eigenvalues, eq.(A.22) of [1], and 
%   the square of the eigenvalues, respectively. LDEX=[min(LD),max(LD)].
%   B2 is (A.26) of [1], and here is of size (number of eigenvalues, number
%   of bits). 
%   POL2 is the number of noise polarizations. NSYMB is the number of bits.
%
%   The function evaluates the saddlepoint by the Newton-Rhapson method
%   along the lines described in [1]. Note that some input parameters are
%   actually not necessary, but they were inserted for speed reasons.
%
%   The accuracy of the routine is taken under control through some
%   parameters at the beginning of the function.
%
%   This function implements the algorithm proposed by E. Forestieri in [1]
%
%   It is recommended to use the equivalent mex file of this
%   function when possible. Try to compile saddle.c by
%
%   mex saddle.c
%
%   If the mex is created successfully, Matlab works with it and discards such
%   function.
%
%   [1] E. Forestieri, "Evaluating the Error Probability in Lightwave 
%       Systems with Chromatic Dispersion, Arbitrary Pulse Shape and Pre-  
%       and Postdetection Filtering,", J. Lightw. Technol., vol. 18, no.11,
%       p. 1493-1503, Nov. 2000.
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


%%%%%%%%%%%%%%%%%
SADTOL = 1e-5;  % tolerance of the newton's method
SADACC = 1e-3;  % accuracy if MAXNX is reached
MAXNX  = 100;    % maximum number of iterations
TOLDEN = 1e-9;	% avoid division by very small numbers
%%%%%%%%%%%%%%%%%

[nno,nsymb] = size(b2);
u0 = zeros(1,nsymb);
nsw = zeros(1,nsymb); % 1: no saddle point. 2: inaccurate saddle point.

for k=1:nsymb

    % Set up the initial guess
    u0(k) = (dx(k)+sgn(k)*sqrt(dx(k)*dx(k)+4.d0*vars))/(2.d0*vars); % (B.11) of [1]
    if u0(k) > 0
        ubi = 1/(dsg2*ldex(2));
        ubl = 0;
        ubr = ubi;
        if (u0(k) >= ubi), u0(k) = 0.99*ubi;end;  % see comment after (B.11) in [1]
    else
        ubl = -Inf;
        ubr = 0;
        if (ldex(1) < 0) 
            ubl = 1/(dsg2*ldex(1));
            if (u0(k) <= ubl), u0(k) = 0.99*ubl;end; % see comment after (B.11) in [1]
        end
    end

    % Search the saddle point with the Newton's method
    nx = 0;
    du = Inf;
    u1 = 0;
    while (abs(du) >= SADTOL) && (nx <= MAXNX) 
        ds2u0 = dsg2*u0(k);
        d1 = -(xi(k)+1/u0(k));
        umbs = 1-ds2u0*ld;
        fnum = b2(:,k)./ld + dsg2*ld.*umbs*pol2;
        fden = umbs.*umbs;
        d1 = d1+sum(fnum./fden);    % first derivative, (B.8) of [1]

        d2 = 1/(u0(k)*u0(k));
        umbs = 1-ds2u0*ld;
        fnum = qsg2*b2(:,k) + qsg4*pol2*ld2.*umbs;
        fden = umbs.*umbs.*umbs;
        d2 = d2 + sum(fnum./fden);  % second derivative, (B.10) of [1]

        if (d2 < eps) 
            if (nsw(k) == 0) 
                nsw(k) = 1;
                return
            end
        end

        nx = nx+1;
        u2 = u1;    % saddle point at cycle nx-2
        u1 = u0(k);    % saddle point at cycle nx-1
        du = -d1/d2;
        u0(k) = u1+du; % saddle point 
        while ((u0(k) <= ubl) || (u0(k) >= ubr))  % Out of bounds?
            du = 0.5*du;	% correct
        u0(k) = u1 + du;
        end	
    end

    if (nx >= MAXNX && (sign(u2-u1) ~= sign(u1-u0(k)))) 
        if (abs((u1-u0(k))/u0(k)) >= SADACC && nsw(k) == 0) 
            nsw(k) = 2;    % low accuracy
        end
    end

    den = u0(k)-2*u1+u2;
    if ((nx >= 2) && (den > TOLDEN)) 
        du = u0(k)-u1;
        u0(k) = u0(k)-du*du/den;  % Aitken's acceleration
    end

end % for k=1:nsymb
