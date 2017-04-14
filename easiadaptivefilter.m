function [Y h1 h2 conv] = easiadaptivefilter(xx, h1, h2, taps, mu, sps )

%EASIADAPTIVEFILTER Source separation filter using EASI algorithm
%   [Y H1 H2] = EASIADAPTIVEFILTER(XX, H1, H2, TAPS, MU, SPS ) applies an 
%   adaptive matrix whose initial state is written in H1 and H2. The 
%   matrix elements are updated with the equivariant adaptive source 
%   separation via independence (EASI) proposed by Cardoso [1]. The 
%   parameters of this algorithm are the convergence parameter MU  and the 
%   number of samples per symbol SPS. TAPS should be always set to 1.
%   Normally there are four elements, let's call them H11, H12, H21, H22 
%   and the relation between the inputs and the outputs of this function is 
%   the following:
%       Y1 = X1 * H11 + X2 * H12
%       Y2 = X1 * H21 + X2 * H22
%   where  X1 and X2 are the two inputs signals and Y1 and Y2 are the 
%   two output signals.
%
%   [1] J.-F. Cardoso, and B. H. Laheld, "Equivariant Adaptive Source 
%   Separation," IEEE Trans. on Signal Processing, vol. 44, no. 12, 
%   Dec 1996
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

N = size( xx, 1 );       % Length of eXtended X = length(X) + taps - 1
ntap = size( h1, 1 );    % Number of Taps
L = N-ntap+1;            % Length of outputs
Y1 = zeros( L , 1 );     % First output vector initialized to zero
Y2 = Y1;                 % Second output vector initialized to zero
%nhtap = floor(ntap / 2);  % (tap number -1) / 2

% This huge for loop is VERY VERY slow in Matlab. It is highly recommended
% to compile the cmaadaptivefilter.c file and to use the mexglx version!
for ktime = 1 : L
    nindex = ktime : ktime + ntap -1;             % time index for inputs
    Y1( ktime  ) =  sum( sum( xx( nindex, : ) .* h1 ) ); % Calculating outputs
    Y2( ktime  ) =  sum( sum( xx( nindex, : ) .* h2 ) );
    Y = [Y1(ktime) Y2(ktime)].';
    E = errorfun( Y , mu);

    h11 = ( 1-mu*E(1,1) ) * h1(:,1) + ( -mu*E(1,2)) * h2(:,1);
    h12 = ( 1-mu*E(1,1) ) * h1(:,2) + ( -mu*E(1,2)) * h2(:,2);
    h21 = ( -mu*E(2,1)  ) * h1(:,1) + (1-mu*E(2,2)) * h2(:,1);
    h22 = ( -mu*E(2,1)  ) * h1(:,2) + (1-mu*E(2,2)) * h2(:,2);        
    
    h1(:,1) = h11;
    h1(:,2) = h12;
    h2(:,1) = h21;
    h2(:,2) = h22;
    
end

Y = [Y1 Y2];

function E = errorfun( Y, lambda)
% Y : m x 1, usually m = 2
% lambda : very small
% E = ( Y*Y' - eye(2) ) ./ (1 + lambda.*abs(Y'*Y)) + ...
%     ( abs(Y).^2.*Y* Y' - Y * (abs(Y).^2.*Y)') ./ (1 + lambda.*abs( Y' * abs(Y)) );
% temp1 = ( Y*Y' - eye(2) ) ./ (1 + lambda.*abs(Y'*Y));
% temp2 = ( abs(Y).^2.*Y* Y' - Y * (abs(Y).^2.*Y)') ./ (1 + lambda.*abs( Y' * abs(Y)) );
a = Y(1);
b = Y(2);
E(1,1) = (abs(a)^2-1)/(1+lambda*(abs(a)^2+abs(b)^2));
E(1,2) = (a*b)/(1+lambda*(abs(a)^2+abs(b)^2))+(a*b*(abs(a)^2-abs(b)^2))/(1+lambda*(a*abs(a)+b*abs(b)));
E(2,1) = (a*b)/(1+lambda*(abs(a)^2+abs(b)^2))+(a*b*(abs(b)^2-abs(a)^2))/(1+lambda*(a*abs(a)+b*abs(b)));
E(2,2) = (abs(b)^2-1)/(1+lambda*(abs(a)^2+abs(b)^2));
