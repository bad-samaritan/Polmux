function y=evaldelay(ftype,bw)

%EVALDELAY Evaluate the group-delay of the filter.
%   Y=EVALDELAY(FTYPE,BW) evaluates the group delay of the filter
%   FTYPE (see MYFILTER for available filters) having 3dB bandwidth 
%   BW.
%
%   Thanks to E. Forestieri for providing the filter expressions.
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


% CONSTANTS

r4p2r2=2.61312592975275;      % =sqrt(2*(2+sqrt(2))); % used in Butterworth 4th order
b1=3.86370330515627315;
Bb=0.3863;

ftype=lower(ftype);
switch ftype
    
    case 'movavg'
% Short term integrator

      y=0;
      
    case 'gauss' 
% Gaussian

      y=0;
      
    case 'gauss_off' 
% Gaussian with offset

      y=0;      
      
    case 'butt2'
% Butterworth 2nd order

      y=1.11*sqrt(2)/(2*pi*bw);
      
    case 'butt4'
% Butterworth 4th order

      y=1.1*r4p2r2/(2*pi*bw);     
    case 'butt6'
% Butterworth 6th order

      y=1.1*b1/(2*pi*bw);
      
    case 'ideal'
% Ideal filter

      y=0;
    
    case 'bessel5'
% Bessel 5th order

      y=Bb/bw;
      
     case 'rc1'
% RC filter 1st order == Butterworth 1st order

      y=1/(2*pi*bw);
      
     case 'rc2'
% RC filter 2nd order

      y=(sqrt(2)-1)/(pi*bw);
      
    case 'supergauss' 
% Super-Gaussian

      y=0;      
    
  otherwise
      error('the filter ftype does not exist.');
end

