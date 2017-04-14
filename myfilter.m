function Hf=myfilter(ftype,f,bw,ord)

%MYFILTER Filter device in the frequency domain.
%   HF=MYFILTER(FTYPE,F,BW,ORD) returns the filter FTYPE evaluated at 
%   frequencies F in the column vector HF. BW is the the 3-dB bandwidth of
%   the filter (with exceptions, Note 2), with respect to the filtered signal 
%   power, i.e., |Hf|^2=0.5 at f=+-bw. 
%   FTYPE can be a string of the following:
%
%   'movavg'    : Short term integrator (moving average) [see Note 2]
%   'gauss'     : Gaussian filter  [see Note 3]
%   'butt2'     : Butterworth 2nd order
%   'butt4'     : Butterworth 4th order
%   'butt6'     : Butterworth 6th order
%   'ideal'     : Ideal filter
%   'bessel5'   : Bessel 5th order
%   'rc1'       : RC1 filter
%   'rc2'       : RC2 filter
%   'supergauss': Super-Gaussian filter of order ORD
%   'gauss_off' : Gaussian filter with offset ORD from the center of the
%                 channel
%
%   Note 1: The bandwidth BW is a lowpass bandwidth -> For optical bandpass
%           filters having 3dB bandwidth Bo, it is BW=Bo/2.
%   Note 2: For the moving average BW is not the 3dB bandwidth,
%           but the first zero of the sinc, i.e. 1/BW is the duration of the
%           moving average. The 3dB bandwidth is 0.443*BW.
%   Note 3: For the Gaussian filter the bandwidth 1/e (B_ue) is related 
%           to the 3dB bandwidth by  B_ue = 1.6986*BW
%   Note 4: if the frequency is normalized, the 3dB bandwidth must be
%           normalized as well to the same value.
%   Note 5: For future new implementations of special filters, use the 
%	    variable ORD for the new filter parameters.
%
%   Thanks to E. Forestieri for the filter expressions.
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
      b2=7.4641016151377546;
      b3=9.1416201726856413;
      b4=b2;
      b5=b1;
      Bb=0.3863;
      d0=945;
      d1=945;
      d2=420;
      d3=105;
      d4=15;

x = f(:)/bw;     % frequency normalized to the bandwidth
ftype = lower(ftype);
switch ftype
    
    case 'movavg'
% Short term integrator

      Hf = sinc(x);
      
    case 'gauss' 
% Gaussian

      Hf = exp(-0.5*log(2).*x.*x);
      
    case 'gauss_off'
% Gaussian with offset

      Hf = exp(-0.5*log(2).*(x-ord/bw).*(x-ord/bw));
      
    case 'butt2'
% Butterworth 2nd order

      Hf = 1./(1-x.*x + i*sqrt(2)*x);
      
    case 'butt4'
% Butterworth 4th order

      x2 = x.*x;
      umx2 = 1-x2;
      Hf = 1./(umx2.*umx2-sqrt(2)*x2 + i*r4p2r2.*x.*umx2);
      
    case 'butt6'
% Butterworth 6th order

      x2 = x.*x;
      x3 = x2.*x;
      x4 = x3.*x;
      x5 = x4.*x;
      x6 = x5.*x;
      Hf = 1./(1.-b2*x2+b4*x4-x6 + i*(b1*x-b3*x3+b5*x5));
      
    case 'ideal'
% Ideal filter

      Hf = (abs(x) <= 1);
    
    case 'bessel5'
% Bessel 5th order

      om  = 2*pi*x*Bb;
      om2 = om.*om;
      om3 = om2.*om;
      om4 = om3.*om;
      om5 = om4.*om;

      pre = d0-d2*om2+d4*om4;
      pim = d1*om-d3*om3+om5;
      Hf = d0./(pre + i*pim);
      
     case 'rc1'
% RC filter 1st order == Butterworth 1st order

      Hf = 1./(1 + i*x);
      
     case 'rc2'
% RC filter 2nd order

      Hf = 1./(1 + i*sqrt(sqrt(2)-1)*x).^2;     % |Hf|^2=0.5 for x=+-1
      
    case 'supergauss'
% Super-Gaussian of order ORD        
      if nargin ~= 4
          error('missing superGauss order');
      end
      
      Hf = exp(-0.5*log(2).*x.^(2*ord));      
    
  otherwise
      error('the filter ftype does not exist.');
end

