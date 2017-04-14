function y=nmod(A,N)

%NMOD N-modulus of an integer.
%   Y=NMOD(A,N) reduces the integer A into 1->N, mod N.
%
%   E.g. N=8. 
%
%   A   ... -2 -1 0 1 2 3 4 5 6 7 8 9 10 ...
%   Y   ...  6  7 8 1 2 3 4 5 6 7 8 1 2  ...
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


y = mod(A-1-N,N)+1;
