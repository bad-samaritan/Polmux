function [Y h1 h2] = cmaadaptivefilter(xx, h1, h2, taps, mu, R, sps)

%CMAADAPTIVEFILTER Polarization demultiplexing filter using CMA algorithm
%   [Y H1 H2] = CMAADAPTIVEFILTER(XX, H1, H2, TAPS, MU, R, SPS) applies a 
%   matrix of adaptive filters whose initial state is written in H1 and H2
%   matrices. The filters coefficients are updated with the constant 
%   modulus algorithm (CMA) proposed by Godard [1]. The parameters of this 
%   algorithm are the number of taps TAPS, the radius R, the convergence 
%   parameter MU and the number of samples per symbol SPS.
%   Normally there are four filters, let's call them H11, H12, H21, H22 and
%   the relation between the inputs and the outputs of this function is the
%   following:
%       Y1 = X1 ** H11 + X2 ** H12
%       Y2 = X1 ** H21 + X2 ** H22
%   where ** is the convolution, X1 and X2 are the two inputs signals
%   and Y1 and Y2 are the two output signals.
%   Usually the filters Hyx have several coefficients (taps).
%   In order to store all the filter a 3-dimensional matrix is needed. I
%   propose the following order for the dimensions: H( tap, y , x ) where H
%   has the following size: 5x2x2 if the filters have 5 taps each.
%   Since Matlab is very slow with operations with such a kind of matrix I
%   preferred to split H in two 2-dimensional matrices: h1 and h2 are the
%   matrix you could virtually obtain by doing: 
%       H1 = squeeze( H( :, 1, : ) )
%       H2 = squeeze( H( :, 2, : ) )
%   So H11 = H1(:,1), H12 = H1(:,2), H21 = H2(:,1), H22 = H2(:,2)
%
%   [1] D. N. Godard, "Self-Recovering Equalization and Carrier Tracking in
%   Two-Dimensional Data Communication Systems," IEEE Trans. Commun., vol. 
%   COM-28, no. 11, Nov. 1980
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
r1=R(1);r2=R(2);
% This huge for loop is VERY VERY slow in Matlab. It is highly recommended
% to compile the cmaadaptivefilter.c file and to use the mexglx version!
for ktime = 1 : L
    nindex = ktime : ktime + ntap -1;             % time index for inputs
    Y1( ktime  ) =  sum( sum( xx( nindex, : ) .* h1 ) ); % Calculating outputs
    Y2( ktime  ) =  sum( sum( xx( nindex, : ) .* h2 ) );
    incr1 = mu.*errorfuncma(Y1(ktime),r1).*conj( xx(nindex,:) ); %Calculating errors
    incr2 = mu.*errorfuncma(Y2(ktime),r2).*conj( xx(nindex,:) );
    h1 = h1 + incr1;                                        % Updating filters
    h2 = h2 + incr2;
end
Y = [Y1 Y2];

function E = errorfuncma( X, M)
E = X.*(M - abs( X ) .^ 2);
