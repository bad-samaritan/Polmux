function create_field(ftype,sigx,sigy,options)

%CREATE_FIELD create the electric field
%   CREATE_FIELD(FTYPE,SIGX) creates the electric field. 
%
%   THIS FUNCTION MUST BE CALLED IN EACH SIMULATION
%
%   SIGX is a matrix [Nfft,Nch] containing on columns the x polarization of
%   the electric fields to be multiplexed togehter. Nch is the overall 
%   number of channels, Nfft the number of FFT points. 
%   CREATE_FIELD creates new fields of the global variable GSTATE, 
%   GSTATE.FIELDX and its copy GSTATE.FIELDX_TX, respectively.
%   GSTATE.FIELDX is the electric field that will be propagated in the 
%   optical system. 
%   
%   FTYPE can be 'sepfields' or 'unique'. With 'sepfields' GSTATE.FIELDX is
%   a copy of SIGX. 'sepfields' is useful for propagation in optical fibers
%   in absence of four wave mixing (FWM) and allows to obtain faster runs.
%   'unique' combines all channels into a unique channel and allows to
%   account for FWM in optical fibers. 'unique' allows therefore more
%   accurate results even if it is slow since requires larger value of
%   GSTATE.NT (see RESET_ALL function) for accounting all channels. With 
%   option 'unique' CREATE_FIELD acts as an ideal multiplexer.
%
%   CREATE_FIELD(FTYPE,SIGX,SIGY) operates on the x (SIGX) and y (SIGY)
%   polarization creating GSTATE.FIELDX, and GSTATE.FIELDY (and their
%   copies, GSTATE.FIELDX_TX and GSTATE.FIELDY_TX, respectively). SIGX and
%   SIGY must have the same size [Nfft,Nch].
%
%   CREATE_FIELD(FTYPE,SIGX,SIGY,OPTIONS) accepts the optional parameter
%   OPTIONS, containing:
%   
%   OPTIONS.DELAY: can be the string 'rand' or a vector of double.
%                  In the first case a uniform distributed random delay
%                  between [0,1] is added to the channels before creating  
%                  the electric field. In the second case the vector is
%                  used as delay. The values are normalized to 
%                  1/GSTATE.SYMBOLRATE, i.e. the symbol time.
%   OPTIONS.POWER: if set to 'average', the power defined in LASERSOURCE is
%                  not the peak power (default) but the average power.
%
%   Note: for hybrid symbol-rate systems, the delay is normalized to the
%   current (last defined) symbol time (1/GSTATE.SYMBOLRATE) and not to 
%   the channel symbol time.
%
%   In absence of y polarization simply use:
%                  CREATE_FIELD(FTYPE,SIGX,[],OPTIONS)
%
%   CREATE_FIELD initializes the global variable GSTATE.DELAY to zero or to
%   the value imposed by OPTIONS.DELAY.
%   CREATE_FIELD initializes the global variable GSTATE.DISP to zero.
%
%   See also: RESET_ALL, LASERSOURCE, ELECTRICSOURCE
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
global CONSTANTS

CLIGHT = CONSTANTS.CLIGHT;

Nfft = GSTATE.NSYMB*GSTATE.NT; % number of FFT points

if exist('options','var')
    checkfields(options,{'delay','power'});
end
[nrx,ncx] = size(sigx);
if nargin >= 3
    if isstruct(sigy)
        error('missing sigy. Set to empty if it does not exist');
    end
    npol = 2; % number of polarizations
    if isempty(sigy) 
        isy = false;        
        npol = 1;
    else
        isy = true;
    end
    if isempty(sigx) 
        error('empty x component');
    end    
    [nry,ncy] = size(sigy);
    if (nry ~= 0) && (ncx ~= ncy || nrx ~= nry)
        error('sigx and sigy must have the same size');
    end
else
    npol = 1;
    isy = false;
    ncy = 0;
end
if (ncx ~= GSTATE.NCH) || (ncy ~= 0 && ncy ~= GSTATE.NCH)
    error('the number of columns of sigx,sigy must be equal to the number of channels');
end

% Set the power to average if necessary 
if exist('options','var') && any(strcmp(fieldnames(options),'power')) && strcmpi(options.power,'average')
    if isy
        avge = mean(abs(sigx).^2+abs(sigy).^2);    % <- 1xNch vector
    else
        avge = mean(abs(sigx).^2);    % <- 1xNch vector
    end
    sigx = sigx .* ( ones(Nfft,1) * ( sqrt(GSTATE.POWER ./ avge )) );
    if isy
        sigy = sigy .* ( ones(Nfft,1) * ( sqrt(GSTATE.POWER ./ avge )) );
    end
    GSTATE.POWER = GSTATE.POWER.*GSTATE.POWER ./ avge;
end

% Now set the delay
if exist('options','var') && any(strcmp(fieldnames(options),'delay'))
    isd = true;
    if strcmp(options.delay,'rand')
        tau1 = rand(npol,GSTATE.NCH);
        tau = round(tau1*GSTATE.NT);
    elseif isnumeric(options.delay)
        [nrd,ncd] = size(options.delay);
        if nrd ~= npol || ncd ~= GSTATE.NCH
            error('the delay must be of size [number of polarizations,number of channels]');
        end
        tau = round(options.delay*GSTATE.NT);
    end
	for kch=1:GSTATE.NCH
        sigx(:,kch) = fastshift(sigx(:,kch),tau(1,kch));
        if isy
            sigy(:,kch) = fastshift(sigy(:,kch),tau(2,kch));
        end
	end	
    GSTATE.DELAY = tau;
else
    GSTATE.DELAY = zeros(npol,GSTATE.NCH);
    isd = false;
end

% Now set the cumulated dispersion
GSTATE.DISP = zeros(npol, GSTATE.NCH);

% Now create the electric fields

if strcmpi(ftype,'sepfields')       % SEPARATE FIELDS
    GSTATE.FIELDX_TX = sigx;
    GSTATE.FIELDX = sigx;
    if isy
        GSTATE.FIELDY_TX = sigy;
        GSTATE.FIELDY = sigy;
    end 
    
elseif strcmpi(ftype,'unique')      % UNIQUE FIELD
    
    lamt = GSTATE.LAMBDA;
    maxl = max(lamt);
    minl = min(lamt);
    fnyqmin = (CLIGHT/minl-CLIGHT/maxl)/GSTATE.SYMBOLRATE;
    if (GSTATE.NT < fnyqmin) && fnyqmin ~= 0
        reply = input('create_field.m: The number of samples per symbol is too small to avoid aliasing. Do you want to continue? Y/N [N]: ', 's');
        if isempty(reply)
            reply = 'N';
        end
        if ~strcmpi(reply,'Y')
            error('number of samples per symbol is too small');
        end
    end

    lamc = 2*maxl*minl/(maxl+minl); %central wavelength: 1/lamc = 0.5(1/maxl+1/minl)
    deltafn = CLIGHT*(1/lamc-1./lamt);   % absolute frequency spacing [GHz]
    minfreq = GSTATE.FN(2)-GSTATE.FN(1);    % minfreq=1/GSTATE.NSYMB
    ndfn = round(deltafn./GSTATE.SYMBOLRATE/minfreq);  % spacing in points
    GSTATE.FIELDX = zeros(Nfft,1);
    zfieldx = fft(sigx);
    if isy, GSTATE.FIELDY = zeros(Nfft,1);end
    for kch=1:GSTATE.NCH    % create the unique field combining the channels       
        GSTATE.FIELDX = GSTATE.FIELDX + fastshift(zfieldx(:,kch),-ndfn(kch));
        if isy
            zfieldy = fft(sigy(:,kch));
            GSTATE.FIELDY = GSTATE.FIELDY + fastshift(zfieldy,-ndfn(kch));
        end
    end
    GSTATE.FIELDX = ifft(GSTATE.FIELDX);
    GSTATE.FIELDX_TX = GSTATE.FIELDX;
    if isy
        GSTATE.FIELDY = ifft(GSTATE.FIELDY);
        GSTATE.FIELDY_TX = GSTATE.FIELDY;
    end
else
    error('ftype must be ''sepfields'' or ''unique''');
end

if GSTATE.PRINT
    
    %%%% PRINT summary

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===           create_field           ===\n');
    fprintf(fid,'========================================\n\n');

    fprintf(fid,'Created GSTATE.FIELDX with option: %s\n',ftype);
    if isy
        fprintf(fid,'Created GSTATE.FIELDY with option: %s\n',ftype);
    end
    fprintf(fid,'\nCreated the copy ..._TX\n');
    if isd && strcmp(options.delay,'rand')
        fprintf(fid,'Added random delay to channels\n\n');
        for kch=1:GSTATE.NCH
            fprintf(fid,'%.3f ',tau1(kch));
        end
        fprintf('\n');
    end
    if isd && isnumeric(options.delay)
        fprintf(fid,'Added user-defined delay to channels\n\n');
        for kch=1:GSTATE.NCH
            fprintf(fid,'%.3f ',options.delay(kch));
        end
        fprintf('\n');

    end    
    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT
    

