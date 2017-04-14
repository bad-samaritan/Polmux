function y=fastshift(x,n)

%FASTSHIFT Fast but simplified circular shift.
%   Y=FASTSHIFT(X,N) shifts the vector or matrix X circularly.
%   If X is a vector then Y=FASTSHIFT(X,N) is equivalent to
%   Y=CIRCSHIFT(X,N). 
%
%   E.g. N=2. 
%
%   X   0 1  2 3 4 5 6 7 8 9 10
%   Y   9 10 0 1 2 3 4 5 6 7  8 
%
%   If X is a matrix of then FASTSHIFT acts on the first dimension
%   (rows) and Y=FASTSHIFT(X,N) is equivalent to Y=CIRCSHIFT(X,[N 0]).
%
%   E.g. N=2. 
%   X = [     1     2     3
%             4     5     6
%             7     8     9 ]
%   Y = [     4     5     6
%             7     8     9
%             1     2     3 ]
%
%   See also CIRCSHIFT, NMOD
%
%   Author: Massimiliano Salsi, 2009
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

if length(n) > 1,error('the shift must be a single element');end
if size(x,1) == 1
    x=x(:);
    spy=true;
else
    spy=false;
end
if n<0
    n=-n;
    y = [x( nmod(n+1,end) :end,:);x(1:mod(n,end),:)];
elseif n>0
    y = [x( nmod( end-n+1,end):end ,:);x( 1:mod(end-n,end),:)];
else
    y=x;
end
if spy
    y=y.';
end
