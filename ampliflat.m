function ampliflat(x,atype,options)

%AMPLIFLAT Ideal Optical amplifier with ASE noise.
%   AMPLIFLAT(X,ATYPE) amplifies the optical field. ATYPE is a string     
%   equal to 'gain' if the amplifier has a flat power gain equal to X [dB].
%   Otherwise, ATYPE can be 'fixpower' if the amplifier takes the gain so
%   as to have an output average power for channel ceil(Nch/2) equal 
%   to X [mW], Nch being the number of channels.
%   This options works only with channels separated (see CREATE_FIELD).
%
%   AMPLIFLAT(X,ATYPE,OPTIONS) has the additional variable OPTIONS to
%   insert the amplified spontaneous emission (ASE) noise.
%   OPTIONS is a structure whose fields can be:
%
%   OPTIONS.f:     [dB] is the optical ASE noise figure, which corresponds 
%                  to a one-sided ASE power spectral density, on two 
%                  polarizations, N0 = F*(Gain-1)*h*nu, with Gain the 
%                  amplifier gain, h the Planck's constant and nu the 
%                  channel central frequency.
%                  Hence, ASE power on a frequency band B is Pase = N0*B.
%   OPTIONS.asepol:If it's 'asex' allows to force to zero the ASE noise 
%                  added to GSTATE.FIELDY, while for 'asey' allows to force
%                  to zero the noise added to GSTATE.FIELDX.
%   OPTIONS.noise: A matrix containing user's defined complex, unit
%                  variance and zero mean Gaussian ASE noise samples. 
%                  OPTIONS.noise must have the same size of 
%                  [GSTATE.FIELDX, GSTATE.FIELDY] and must be read in that 
%                  way.
%
%   If the amplifier does not generate ASE noise, don't set OPTIONS.
%   
%   Note: AMPLIFLAT assumes the same gain for both polarizations. 
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

global CONSTANTS;  % physical constants: a structure defined in reset_all.m
global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

HPLANCK = CONSTANTS.HPLANCK;      % Planck's constant [J*s]
CLIGHT = CONSTANTS.CLIGHT;      % speed of light [m/s]

[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 

atype = lower(atype);
switch atype
    case 'gain'
        gain = 10^(x*0.1);
    case 'fixpower'                 % gain -> out power = x for ch.1
        if nfc == GSTATE.NCH
            midch = ceil(nfc/2);
            avge = avg_power(midch,'abs'); % average energy per bit ch.midch
            gain = x/avge;    
        else
            error(['''fixpower'' works',...
                ' only for channels separated']);
        end
                                                    
    otherwise
        error('wrong string atype');
end
GSTATE.FIELDX = GSTATE.FIELDX*sqrt(gain);
isy = ~isempty(GSTATE.FIELDY);
if isy
    GSTATE.FIELDY = GSTATE.FIELDY*sqrt(gain);
end

%%% Add the noise

% noy1 = 2;   % 2: no ASE. 1: ASE on X and Y pol. 0: ASE on X pol.

if (nargin == 3)
    fields = fieldnames(options);
    noy1 = 1;
    if (max(strcmp(fields,'f'))  &&  ~isinf(options.f))
        Flin = 10^(options.f*0.1);
        if nfc==1
            maxl=max(GSTATE.LAMBDA);
            minl=min(GSTATE.LAMBDA);
            lamc = 2*maxl*minl/(maxl+minl);     % central wavelength
            sigma = sqrt(Flin/4*HPLANCK*CLIGHT./lamc.*(gain-1)*...
            GSTATE.NT.*GSTATE.SYMBOLRATE*1e21);   % sqrt(mW). sigma^2 is the variance
        else
            sigma = sqrt(Flin/4*HPLANCK*CLIGHT./GSTATE.LAMBDA.*(gain-1)*...
            GSTATE.NT.*GSTATE.SYMBOLRATE*1e21);   % sqrt(mW). sigma^2 is the variance
        end
    else
        sigma=0;
    end
    if sigma
        if max(strcmp(fields,'onepol'))
            if strcmpi(options.onepol,'asex')   % noise only on X polarization
                noy1 = 0; % only x
                asepol = [true false];
            elseif strcmpi(options.onepol,'asey')
                noy1 = -1; % only y
                asepol = [false true];
            else
                error('ONEPOL, if exists, must be ''asex'' or ''asey''');
            end
        else
            asepol = [true true];
        end

        sigma_matrix = ones(nfr,1) * sigma;    % size nfr x Nch

        if max(strcmpi(fields,'noise'))  % No need to generate noise
            if asepol(1)
                noiseX = sigma_matrix .* options.noise(:,1:nfc);
            end
            if asepol(2) % here the first nfc columns exist and are zeros
                noiseY = sigma_matrix .* options.noise(:,nfc+1:end);
            end
        else
            if asepol(1)
                noiseX = sigma_matrix .* (randn(nfr,nfc)+i*randn(nfr,nfc));
            end
            if asepol(2)
                noiseY = sigma_matrix .* (randn(nfr,nfc)+i*randn(nfr,nfc));
            end
        end
        if asepol(1)
            GSTATE.FIELDX = GSTATE.FIELDX + noiseX;
        end
        if asepol(2)
            if isy
                GSTATE.FIELDY = GSTATE.FIELDY + noiseY;
            else
                GSTATE.FIELDY = noiseY;
                GSTATE.DELAY(2,:) = zeros(1,GSTATE.NCH);
            end
        end
        
    end
else
    noy1 = 2;
end


if GSTATE.PRINT
    %%%% PRINT summary

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===            ampliflat             ===\n');
    fprintf(fid,'========================================\n\n');
    switch atype
	case 'gain'
            fprintf(fid,'Amplifier Gain: %.2f  [dB]  (linear: %.2f)\n',...
        	x,gain);
	case 'fixpower'
            fprintf(fid,'Amplifier with user defined output power\n');
            fprintf(fid,'correspondent to a gain %.2f [dB] (linear: %.2f)\n',...
        	10*log10(gain),gain);
    end
    if noy1  == 2
        addstr = 'without';
        fprintf(fid,'\nNO ASE noise');
    else
        addstr = 'with';
        if isempty(options.f), options.f = -Inf;end;
        if noy1 == 1
            fprintf(fid,'\nAdded ASE noise on both Signal polarizations with:\n');
        elseif noy1 == 0
            fprintf(fid,'\nAdded ASE noise ONLY on X polarization with:\n');
        else
            fprintf(fid,'\nAdded ASE noise ONLY on Y polarization with:\n');            
        end
        fprintf(fid,'Noise figure: %.2f [dB]\n',options.f);
    end
    if nfc == GSTATE.NCH
        fprintf(fid,['\nOutput Average Power, ',addstr,' noise:\n\n']);
        for kch=1:GSTATE.NCH
            avgp = avg_power(kch,'abs');
            fprintf(fid,'ch. #%d: %8.4f [mW] (%6.2f [dBm])\n',kch,...
                avgp,10*log10(avgp));
        end
    end

    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
end
    
        
       
