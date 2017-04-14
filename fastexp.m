function y=fastexp(x)
%FASTEXP Calculates exp( i * x ) quickly
%   Y=FASTEXP(X) 
%   It's a fast computation of the expression:
%       Y = EXP( i * X )
%   If you have compiled the file fastexp.c with mex then Matlab will run
%   the still quicker (or quickest) fastexp.mexglx
%
%   Author: Massimiliano Salsi, 2009
%	Author: Paolo Serena, 2009
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
y=complex(cos(x),sin(x));
