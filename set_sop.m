function [Ex,Ey,matout]=set_sop(Ex,Ey,ang1,ang2,angtype,rtype)

% SET_SOP sets the average State Of Polarization of a signal.
%   [EX,EY]=SET_SOP(EX,EY,ANG1,ANG2) sets the State Of Polarization  
%   (SOP) of the electric field having x component EX and y component EY. 
%   EX and EY are column vectors (see MZ_MODULATOR). ANG1, ANG2 are the
%   azimuth and ellipticity [rad] of the output SOP.
%
%   [EX,EY]=SET_SOP(EX,EY,ANG1,ANG2,ANGTYPE) specifies ANG1 and ANG2 
%   according to ANGTYPE. Available options are:
%
%   ANGTYPE = 'aziell'; user specifies azimuth and ellipticity angles.
%       Latitude and longitude on the Poincaré sphere equal 2*ANG2 and
%       2*ANG1, respectively.
%   NOTE: azimuth, normally in [-pi/2;pi/2], is interpreted modulo pi;
%         ellipticity, normally in [-pi/4;pi/4], is interpreted modulo pi/2;
%
%   ANGTYPE = 'aarphd'; user specifies the "absolute amplitude ratio"
%       atan(abs(EY/EX)) and the phase difference between field components
%       arg(EY*conj(EX)).
%   NOTE: absolute amplitude ratio, normally in [0;pi/2], is interpreted 
%         modulo pi/2;
%         phase difference, normally in [-pi;pi], is interpreted modulo
%         2*pi.
%
%   [EX,EY,MAT]=SET_SOP(EX,EY,ANG1,ANG2) also returns in MAT the unitary
%   matrix that rotates the SOP. MAT can be used by INVERSE_PMD.
%
%   [EX,EY]=SET_SOP(EX,EY,ANG1,ANG2,RTYPE) with RTYPE='mean' change the 
%   default behavior by setting the output SOP with ANG1 and ANG2 rotations 
%   relatively to a reference system aligned with the input average SOP. 
%   The degree of polarization (DOP) is preserved.
%
%   Example: an ideal polarization division multiplexed signal (PDM) signal  
%       before SET_SOP lies in the plane (S2,S3). 
%       Calling [EX,EY]=SET_SOP(EX,EY,pi/2,pi/4) let it lie in the plane 
%       (S1,S2).
%
%   See also: DOP_METER, MZ_MODULATOR, INVERSE_PMD
%
%   Author: Armando Vannucci, 2009
%   Contributed: Paolo Serena, 2009
%                Marco Bertolini, 2009
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



if ~exist('Ey','var'), 
    error('trying to set SOP of scalar optical field... absurd'); 
end
if ~exist('angtype','var'), angtype='aziell';end   % default

% get the average SOP
% Evaluate Stokes coordinates, for every sample of k-th channel
if exist('rtype','var') 
    if strcmp(rtype,'mean')
        S0 = mean(abs(Ex).^2 + abs(Ey).^2);
        S1 = mean(abs(Ex).^2 - abs(Ey).^2);
        S2 = mean( 2.*real( Ex.*conj(Ey)));
        S3 = mean(-2.*imag( Ex.*conj(Ey))); % ref system: average SOP
    else
        error('RTYPE must be ''mean''');
    end
else
    S0 = 1; S1 = 1; S2 = 0; S3 = 0;     % ref system: S1
end
    % get azimuth and ellipticity of the average SOP
    azi= sign(S2).*acos(S1./sqrt(S0.^2-S3.^2))/2 + (S2==0).*(1-sign(S1)).*pi/4;
    ell= asin(S3./S0)/2;

% matrix that brings the average SOP onto Stokes S1 axis
matS1 = ...
   [ cos(azi)*cos(ell)+i*sin(azi)*sin(ell), ...
     sin(azi)*cos(ell)-i*cos(azi)*sin(ell); ...
    -sin(azi)*cos(ell)-i*cos(azi)*sin(ell), ...
     cos(azi)*cos(ell)-i*sin(azi)*sin(ell) ];
switch angtype
    case 'aziell'
        % matrix that brings Stokes S1 axis onto the selected SOP
        matSOP = ...
           [ cos(ang1)*cos(ang2)-i*sin(ang1)*sin(ang2), ...
            -sin(ang1)*cos(ang2)+i*cos(ang1)*sin(ang2); ...
             sin(ang1)*cos(ang2)+i*cos(ang1)*sin(ang2), ...
             cos(ang1)*cos(ang2)+i*sin(ang1)*sin(ang2) ];
    case 'aarphd'
        % matrix that brings Stokes S1 axis onto the selected SOP
        matSOP = ...
           [ cos(ang1)*cos(ang2/2)-i*cos(ang1)*sin(ang2/2), ...
            -sin(ang1)*cos(ang2/2)+i*sin(ang1)*sin(ang2/2); ...
             sin(ang1)*cos(ang2/2)+i*sin(ang1)*sin(ang2/2), ...
             cos(ang1)*cos(ang2/2)+i*cos(ang1)*sin(ang2/2) ];
    otherwise
        error('angles type must be either ''aziell'' or ''aarphd'' ');
end

% overall matrix to apply
mat= matSOP*matS1;
% to bring the field onto the selected SOP
dummyX= Ex; dummyY= Ey;   % make a local copy
Ex = mat(1,1)*dummyX + mat(1,2)*dummyY;
Ey = mat(2,1)*dummyX + mat(2,2)*dummyY;
if nargout == 3, 
    matout = mat; 
end
