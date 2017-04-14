function E=lasersource(Ptx,lam,spac,options)

%LASERSOURCE Multichannel laser transmitter
%   E = LASERSOURCE(PTX,LAM,SPAC,OPTIONS) creates the WDM optical field 
%   whose channels are saved into the columns of the matrix E. 
%
%	PTX contains the channel's peak power. PTX can be vector [1,Nch], being 
%	Nch the number of channels, or a scalar.
%   In the last case, the same value is used for all channels. 
%	PTX is saved by this function into the global variable GSTATE.POWER.
%
%   LAM are the wavelengths [nm] of the channels. In this function they are
%   associated with the global variable GSTATE.LAMBDA, here and forever.
%   LAM can be a vector [1,Nch] or a scalar: if it is a scalar,
%   the additional parameter SPAC is required, which indicates the spacing
%   between channels [nm], while LAM is assumed as the central wavelength.
%
%   E.g.  Nch=4, lam=1550, spac=0.8 -> lambda=[1548.8 1549.6 1550.4 1551.2]
%
%   OPTIONS is a structure whose (optional) fields can be:
%
%      OPTIONS.single: if exists and it is true, LASERSOURCE generates a
%           single laser line and thus, PTX, LAM and LINEWIDTH, N0, if  
%           present, must be scalar. In this case SPAC must still be  
%           specified, but its value is neglected and thus could be safely 
%           set to 0.
%      OPTIONS.linewidth: represents the 3 dB width of the laser line, 
%           normalized to the symbolrate. LINEWIDTH can be a scalar or a 
%           vector whose length equals the number of channels
%      OPTIONS.n0: represents the one-sided spectral density of a gaussian 
%           complex noise added to the laser, in dB. This way it's possible  
%           to set the desired OSNR of the laser.
%      OPTIONS.anoise: the matrix with amplitude noise samples (real and 
%           imag).
%      OPTIONS.pnoise: the matrix with phase noise samples
%
%   The matrices representing amplitude and phase noise are ignored if
%   OPTIONS.linewidth and OPTIONS.n0 are specified.
%
%   Note: GSTATE.POWER may be modified by CREATE_FIELD if the user 
%       indicates that the power is the average power. See CREATE_FIELD.
%
%   See also MZ_MODULATOR, ELECTRICSOURCE, PATTERN, CREATE_FIELD
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


global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m
% NOTE: In this function the field GSTATE.LAMBDA and GSTATE.POWER.
% will be initialized

nfft = GSTATE.NSYMB*GSTATE.NT;
nch = GSTATE.NCH;

if nargin == 4      % Checking options
    optfields = fieldnames(options);
    checkfields(options,{'n0','linewidth','anoise','pnoise','single'});
    if max(strcmpi(optfields,'n0'));
        N0 = options.n0;
    else
        N0 = -inf;
    end

    if max(strcmpi(optfields,'linewidth'));
        linewidth = options.linewidth;
    else
        linewidth=0;
    end

    if max(strcmpi(optfields,'anoise'))
        anoise = options.anoise;
    end

    if max(strcmpi(optfields,'pnoise'))
        pnoise = options.pnoise;
    end

    if max(strcmpi(optfields,'single'))
        single = options.single;
    else
        single = false;
    end
else
    single = false;
    linewidth = 0; 
    N0 = -Inf;
end     % end of options checks

% check that the length of the amplitude noise
if length(N0) ~= 1
    error('N0 length must be 1');
end

if single
    if length(linewidth) ~= 1
        error(['you selected the ''single'' ',...
            'option, linewidth length must be 1']);
    end
    if length(Ptx) ~= 1
        error(['you selected the ''single'' ',...
            'option, Ptx length must be 1']);
    end 
    if length(lam) ~= 1
        error(['you selected the ''single'' ',...
            'option, lam length must be 1']);
    end        
    
    if sum(GSTATE.LAMBDA==lam)
        reply = input(['You are putting two channels on the same wavelength. ',...
            'Do you want to continue? Y/N [N]: '], 's');
        if isempty(reply)
            reply = 'N';
        end
        if ~strcmpi(reply,'Y')
            error('putting two channels on the same wavelength');
        end
    end
    GSTATE.LAMBDA = [GSTATE.LAMBDA lam];
    GSTATE.POWER = [GSTATE.POWER Ptx];
    Pin = Ptx;
    nch = 1;
else

    % check that the length of the linewidth
    if (length(linewidth) ==1 && nch > 1)
        linewidth=linewidth*ones(1,nch);
    elseif length(linewidth)~=nch
        error('linewidth length must be 1 or Nch');
    end

    Lp = length(Ptx);  % init power
    if Lp == 1
        Pin = Ptx*ones(1,nch);
    elseif Lp ~= nch
        error('wrong length for PTX (must be 1 or # channels)');
    else
        Pin = Ptx;
    end
    La = length(lam);   % init wavelengths
    if (La == 1) && (nch > 1)
        if ~exist('spac','var')
            error('missing the channel-spacing SPAC');
        end
        lamt = zeros(1,nch);
        for chan = 1:nch
            lamt(chan)=lam+spac*(chan-(nch+1)/2);
        end
    elseif La ~= nch;
        error('wrong length for LAM (must be 1 or # channels)');
    else
        lamt = lam;
    end
    
    maxl = max(lamt);
    minl = min(lamt);
    GSTATE.LAMBDA = lamt;   % global and unique definition of wavelengths
    lamc = 2*maxl*minl/(maxl+minl); %central wavelength: 1/lamc = 0.5(1/maxl+1/minl)

    GSTATE.POWER = Pin;

end % end if single

E   = ones(nfft,1)*sqrt(Pin);

% Add phase noise
if linewidth
    freq_noise  = (ones(nfft,1) * sqrt(2*pi*linewidth./GSTATE.NT)) .* randn( nfft, nch);
    freq_noise(1) = 0;
    phase_noise = cumsum(freq_noise,1);
    % Brownian bridge
    for nnoise=1:length(phase_noise)
        phase_noise(nnoise)=phase_noise(nnoise)-(nnoise-1)/...
            (length(phase_noise)-1)*phase_noise(end);
    end
    E = E .* fastexp(phase_noise);
elseif exist('pnoise','var')
    E = E .* fastexp(pnoise);
end

% Add gaussian complex white noise (ASE)
if ~isinf(N0)
    N0_lin=10^(.1*N0);
    ampl_noise = ( ones(1,nch) * sqrt(N0_lin) ) .* ...
        complex(randn(nfft,nch),randn(nfft,nch)) ./ sqrt(nch-nfc+1);
    E = E + ampl_noise;
elseif exist('anoise','var')
    E = E + anoise;
end

if GSTATE.PRINT

    %%%% PRINT summary
    pstr = 'peak';
    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===            lasersource           ===\n');
    fprintf(fid,'========================================\n\n');
    if Lp == 1
        fprintf(fid,'%s power: %.4f [mW] (%.2f [dBm])\n',pstr,Pin(1),...
            10*log10(Pin(1)));
        if (nch > 1) && ~single, fprintf(fid,'equal for all channels\n');end;
    else
        for kch=1:nch
            fprintf(fid,'%s power ch. #%d: %.4f [mW] (%.2f [dBm])\n',...
                pstr,kch,Pin(kch),10*log10(Pin(kch)));
        end
    end
    if La == 1
        fprintf(fid,'\nCentral wavelength: %.2f [nm]\n',lamc);
        if (nch > 1) && ~single, fprintf(fid,'Channel spacing: %10.4f [nm]',spac);end;
    else
        for kch=1:nch
            fprintf(fid,'\nWavelength ch. #%d: %.2f [nm]',kch,...
                GSTATE.LAMBDA(kch));
        end
    end
    fprintf(fid,['\n\nWavelengths inserted into the global variable ',...
        'GSTATE.LAMBDA\n\n']);

    if length(linewidth) == 1
        fprintf(fid,'Laser linewidth: %.6f [a. u.] \n',linewidth);
        if (nch > 1) && ~single
            fprintf(fid,'equal for all channels\n');
        end
    else
        for kch=1:nch
            fprintf(fid,'#%d Laser linewidth: %.6f [a. u.] \n',kch,linewidth);
        end
    end
    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);

end % end IF GSTATE.PRINT




