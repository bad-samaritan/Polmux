function pol_scrambler(type,coh_timeR,theta,epsilon,Delphi)
% POL_SCRAMBLER rotates the SOP of signal samples on the Poincare sphere.
%   POL_SCRAMBLER(TYPE, COH_TIMER, THETA, EPSILON,DELPHI)  rotates the signal State Of
%   Polarization (SOP) on the Poincaré sphere. 
%   If TYPE = 'rand' then a random rotation is used, so that the signal SOP
%   is uniformly scattered on the Poincaré sphere, as is typical of "long"
%   SMF fibers. The rotation parameters are changed randomly every
%   COH_TIMER symbol periods.
%   For the method used to randomize rotations, see eqs.(13-15) and
%   following observation in [1].
%   COH_TIMER  is the "coherence time" of the scrambler, normalized to the
%   symbol interval (i.e., multiplied by the symbol rate "R"). Use
%   COH_TIMER= 1/NT to scramble each sample independently and produce a DOP
%   that approaches zero; COH_TIMER= NSYMB to scramble all samples equally
%   and preserve the original DOP (NSYMB=# of transmitted symbols; NT=# of
%   samples per symbol).
%
%   If TYPE = 'fixed' then THETA, EPSILON and DELPHI should be provided:
%   the signal SOP is then rotated by a fixed  amount DELPHI [rad.] (the Jones matrix "retardation")
%   about a fixed rotation axis (the Jones matrix "eigenmode") with azimuth THETA
%   and ellipticity EPSILON. 
%   COH_TIMER is ignored in this case.
%
%   [1] A.Vannucci,A.Bononi, "Statistical Characterization of the Jones Matrix
%   of Long Fibers Affected by Polarization Mode Dispersion (PMD)", IEEE JLT
%   may2002, pp.811-821.
%
%   See also: SET_SOP, DOP_METER
%
%   Author: Armando Vannucci, 2009
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

global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

if strcmpi(type,'rand')
    if nargin>2
        warning('function pol_scrambler: additional parameters will be ignored in random mode.');
    end
    % number of samples in one scrambling period (step) and number of scrambling steps
    n_samp_step= max(1,GSTATE.NT*coh_timeR);    % at least, 1 sample
    if (coh_timeR>GSTATE.NSYMB), 
        n_samp_step= GSTATE.NSYMB*GSTATE.NT;    % at most, all samples are equally scrambled
    end;
    n_scram_steps= ceil(GSTATE.NSYMB*GSTATE.NT/n_samp_step);
    % azimuth and ellipticity of the fiber eigenmode, which is uniform on
    % Poincaré sphere: pdfs taken from [1]
    theta=rand(n_scram_steps,1)*pi - 0.5*pi;     % Uniform azimuth in [-pi/2;pi/2]
    epsilon= 0.5*asin(rand(n_scram_steps,1)*2-1);    % so that f_{epsilon}(e)=cos(2*e) in [-pi/4<e<pi/4]
    % retardation Delphi with pdf [1] f_{Delphi}(p)=(1-cos(p))/(2*pi) in [0<p<2*pi]
    uu= rand(n_scram_steps,1);
    % We apply the "percentile transformation method" to a uniform RV uu;
    % The CDF F_{Delphi}(p)=(p-sin(p))/(2*pi) of the retardation is not analytically
    % invertible, hence we use Newton-Raphson's method to find
    % p=F^{-1}_{Delphi}(uu)
    Delphi= pi*ones(n_scram_steps,1);  % initial guess
    while (max(abs(Delphi-sin(Delphi)-2*pi*uu))>1e-3),   % check convergence
        Delphi= Delphi-(Delphi-sin(Delphi)-2*pi*uu)./(1-cos(Delphi));    % update by linear interpolation
    end;
elseif strcmpi(type,'fixed')
    if nargin<5
        error('function pol_scrambler: theta, epsilon and DeltaPhi parameters are required in fixed mode.');
    end
    n_scram_steps= 1;   % all samples equally scrambled with provided parameters
    n_samp_step= GSTATE.NSYMB*GSTATE.NT;
else
    error('Unknown mode.');
end

% Pauli matrices
sig0  = eye(2);
sig1  = [1 0; 0 -1];
sig2  = [0 1; 1  0];
sig3i = [0 1;-1  0]; % =i*sig3 = i*[0 -i;i 0]

for i_step= 1:n_scram_steps,    % samples scrambled in block of length equal to the scrambler coherence time

    % Jones matrix implementing a rotation of Delphi [rad.] around a Stokes
    % axis with azimuth theta and ellipticity epsilon
    matR= complex(cos(Delphi(i_step)/2)*sig0-sin(Delphi(i_step)/2)*sin(2*epsilon(i_step))*sig3i , -sin(Delphi(i_step)/2)*cos(2*theta(i_step))*cos(2*epsilon(i_step))*sig1-sin(Delphi(i_step)/2)*sin(2*theta(i_step))*cos(2*epsilon(i_step))*sig2);

    if (abs(det(matR)-1) > 1e-9),
       error('function pol_scrambler: scrambling Jones matrix with det[]<>1');
    end;    % if (abs(det(matR)-1) > 1e-9),

    % signal samples to scramble
    first_sample= (i_step-1)*n_samp_step+1; last_sample= min(first_sample+n_samp_step-1,GSTATE.NSYMB*GSTATE.NT);
    % Apply scrambling Jones matrix
    dummyX= GSTATE.FIELDX(first_sample:last_sample); dummyY= GSTATE.FIELDY(first_sample:last_sample);   % make a local copy
    GSTATE.FIELDX(first_sample:last_sample) = matR(1,1)*dummyX + matR(1,2)*dummyY;
    GSTATE.FIELDY(first_sample:last_sample) = matR(2,1)*dummyX + matR(2,2)*dummyY;

end;    % for i_step= 1:n_scram_steps,

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRINT SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GSTATE.PRINT
    outfile = [GSTATE.DIR,'/simul_out'];
    fid = fopen(outfile,'a');
    fprintf(fid,'========================================\n');
    fprintf(fid,'===           pol_scrambler          ===\n');
    fprintf(fid,'========================================\n');
    if (type=='rand')
        fprintf(fid,'SOP rotation with Randomly generated parameters \n');
        fprintf(fid,'coherence time of %d symbols implied %d scrambling steps \n',coh_timeR,n_scram_steps);
        fprintf(fid,'SOP rotated by DPhi [deg] around an axis with azimuth theta [deg] and ellipticity epsilon [deg] \n');
%         fprintf(fid,'DPhi    theta      epsilon \n');
%         fprintf(fid,'%5.0f %6.0f %8.0f \n',[180*Delphi'/pi;180*theta'/pi;180*epsilon'/pi]);
        fprintf(fid,'\n');
    else
        fprintf(fid,'SOP rotation with fixed parameters (provided by user) \n');
        fprintf(fid,'SOP rotated by %3.0f [deg] around an axis with azimuth %2.0f [deg] and ellipticity %2.0f [deg] \n\n',180*Delphi/pi,180*theta/pi,180*epsilon/pi);
    end;
    fclose(fid);       
end % end IF GSTATE.PRINT

return; % function polscrambler()
