function q=ber2q(ber)
%Q=BER2Q(BER) Convert the bit-error rate in Q-factor
%	Q=BER2Q(BER) Converts the bit-error rate in Q-factor [dB]
%	using the formula:
%
%	Q=20*log10( sqrt(2) * erfcinv(2*ber));
%
%	As a reference (BER -> Q [dB}): 
%						1e-3 -> 9.80		8e-4 ~= 10
%						1e-5 -> 12.60		2e-4 ~= 11
%						1e-9 -> 15.56		
%
%   Author: Massimiliano Salsi, 2009
%   Comments: Armando Vannucci, 2009
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

q=20*log10( sqrt(2) * erfcinv(2*ber));

% NOTES ON THE CALCULATION
% Based on the "gaussian approximation", the statistics of the received photocurrent, conditioned on the transmission of a "0" or a "1" bit,
% is a gaussian function with mean values m0 or m1 and std. deviations s0 or s1, respectively.
% Under such assumption, the BER can be expressed as
% BER= Q(q)
% where q is the so-called "Q factor", obtained by setting the decision threshold to the crosspoint of the gaussian pdfs describing 
% the photocurrent distribution conditioned on the transmission of a "0" or a "1": such threshold, which is the optimal choice in 
% the "minimax sense", is th=(m0*s1+m1*s0)/(s1+s0).
% Hence, q=(th-m0)/s0=(m1-th)/s1, i.e., the mean-threshold distance divided by the std. deviation, is also equal to q=(m1-m0)/(s1+s0).
% The latter expression of q resembles a "signal-to-noise ratio", (m1-m0) being the swing of received bit amplitudes 
% and (s1+s0) being the rms amplitude of noise.
%
% Q(x) is defined as the area, from abscissa x to infinity, under a normalized Gaussian function, i.e., with zero mean and unit variance.
% This is calculated resorting to the "complementary error function" erfc() in the math library, which is instead twice the area,
% from abscissa x to infinity, under a gaussian function with zero mean and std. deviation equal to 0.5.
% It is thus simple to find the correspondence between Q() function and efrc():  Q(q)= 0.5*erfc(q/sqrt(2)).
% This is the reason why, in linear scale, BER=Q(q) yields q= sqrt(2)*erfcinv(2*BER).
