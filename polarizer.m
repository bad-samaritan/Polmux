function polarizer(ang1,ang2,angtype,T)

%POLARIZER Linear optical polarizer
%   POLARIZER(ANG1,ANG2) polarizes the electric field contained in
%   GSTATE (see RESET_ALL) into the direction having azimuth ANG1 and
%   ellipticity ANG2.
%
%   POLARIZER(ANG1,ANG2,ANGTYPE) has the additional ANGTYPE selecting the 
%   type of angular coordinates used to specify the position of the SOP on 
%   the Poincaré sphere. Available options are:
%
%   ANGTYPE = 'aziell'; user specifies azimuth and ellipticity angles.
%       Latitude and longitude on the Poincaré sphere equal 2*ANG2 and
%       2*ANG1, respectively.
%   NOTE: azimuth, normally in [-pi/2;pi/2], is interpreted modulo pi;
%         ellipticity, normally in [-pi/4;pi/4], is interpreted modulo pi/2;
%
%   ANGTYPE = 'aarphd'; user specifies the "absolute amplitude ratio"
%       atan(abs(GSTATE.FIELDY/GSTATE.FIELDX)) and the phase difference 
%	between field components arg(GSTATE.FIELDY*conj(GSTATE.FIELDX)).
%   NOTE: absolute amplitude ratio, normally in [0;pi/2], is interpreted 
%         modulo pi/2;
%         phase difference, normally in [-pi;pi], is interpreted modulo
%         2*pi.
%
%   POLARIZER(ANG1,ANG2,ANGTYPE,T) indicates in T = [TM, Tm] the power
%   transmittance of the polarizer, being TM the maximum power flow and Tm
%   the minimum power flow along the polarizer's axes. TM and Tm must be
%   within 0 and 1.
%
%   See also: POL_SCRAMBLER, SET_SOP, DOP_METER, RESET_ALL
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


global GSTATE

if nargin == 2 % default
    angtype = 'aziell';
end
if nargin == 4
    if length(T) ~= 2, error('Vector T must have two elements');end
    if (T(1) > 1) || (T(1) < 0) 
        error('Transmittance T(1) must be in the range 0 <= T(1) <= 1');
    end
    if (T(2) > 1) || (T(2) < 0)
        error('Transmittance T(2) must be in the range 0 <= T(2) <= 1');
    end  
    TM = T(1); Tm = T(2);
else
    TM = 1;
    Tm = 0;
end
switch angtype
    case 'aziell'
        % matrix that brings Stokes S1 axis onto the selected SOP
        matrot = ...
           [ cos(ang1)*cos(ang2)-i*sin(ang1)*sin(ang2), ...
            -sin(ang1)*cos(ang2)+i*cos(ang1)*sin(ang2); ...
             sin(ang1)*cos(ang2)+i*cos(ang1)*sin(ang2), ...
             cos(ang1)*cos(ang2)+i*sin(ang1)*sin(ang2) ];
    case 'aarphd'
        % matrix that brings Stokes S1 axis onto the selected SOP
        matrot = ...
           [ cos(ang1)*cos(ang2/2)-i*cos(ang1)*sin(ang2/2), ...
            -sin(ang1)*cos(ang2/2)+i*sin(ang1)*sin(ang2/2); ...
             sin(ang1)*cos(ang2/2)+i*sin(ang1)*sin(ang2/2), ...
             cos(ang1)*cos(ang2/2)+i*cos(ang1)*sin(ang2/2) ];
    otherwise
        error('angles type must be either ''aziell'' or ''aarphd'' ');
end
matdiag = [sqrt(TM) , 0; 0 , sqrt(Tm)];
mat = matrot*matdiag*matrot';
dummyX= GSTATE.FIELDX; dummyY= GSTATE.FIELDY;   % make a local copy
GSTATE.FIELDX = mat(1,1)*dummyX + mat(1,2)*dummyY;
GSTATE.FIELDY = mat(2,1)*dummyX + mat(2,2)*dummyY;
