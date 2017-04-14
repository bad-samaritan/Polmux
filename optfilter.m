function optfilter(ich,ftype,bw,ord)

%OPTFILTER Optical filter.
%   OPTFILTER(ICH,FTYPE,BW) filters the optical field GSTATE.FIELDX and
%   GSTATE.FIELDY with the filter FTYPE (see MYFILTER) having 3-dB
%   bandwidth BW (normalized to the symbolrate), with central frequency 
%   centered on the ICH channel.
%   The resulting filtered field is again associated to GSTATE.FIELDX 
%   and GSTATE.FIELDY, respectively, hence the the field is overwritten. 
%
%   In the case of separated fields (see RESET_ALL) only the ICH column 
%   of GSTATE.FIELDX and GSTATE.FIELDY is filtered, while the other columns 
%   remain unchanged. The global variable GSTATE.DELAY is updated.
%
%   OPTFILTER(ICH,FTYPE,BW,ORD) use the additional parameter ORD for special
%   filter. E.g. ORD is the supergauss order for the supergaussian filter 
%   (see MYFILTER)
%
%   See Also MYFILTER
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

global CONSTANTS;  % CONSTANTS is a global structure variable.
CLIGHT = CONSTANTS.CLIGHT;      % speed of light in vacuum [m/s]

global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

if (ich < 1) || (ich > GSTATE.NCH)
    error('wrong channel number in ICH');
end
if nargin < 4
    ord = 0;     % not using supergaussian optical filter
end
[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
                                    % for all channels.
isy = ~isempty(GSTATE.FIELDY);           

%%%% FILTER THE CHANNEL

Hf = myfilter(ftype,GSTATE.FN,bw*0.5,ord); % Remember that the in the lowpass
%   equivalent domain, the 3 dB bandwidth goes from -bw/2 to + bw/2
if (nfc == GSTATE.NCH) 
    tpower = GSTATE.POWER(ich);
    Eb_before = avg_power(ich);
    GSTATE.FIELDX(:,ich) = ifft(fft(GSTATE.FIELDX(:,ich)) .* Hf);
    if isy
        GSTATE.FIELDY(:,ich) = ifft(fft(GSTATE.FIELDY(:,ich)) .* Hf);
    end
    Eb_after = avg_power(ich);    
else    
    Nfft = GSTATE.NSYMB*GSTATE.NT;
    minfreq = GSTATE.FN(2)-GSTATE.FN(1);  
    maxl=max(GSTATE.LAMBDA);
    minl=min(GSTATE.LAMBDA);
    lamc = 2*maxl*minl/(maxl+minl); %central wavelength
    deltafn = CLIGHT*(1/lamc-1./GSTATE.LAMBDA(ich)); % frequency spacing
    ndfn = round(deltafn./GSTATE.SYMBOLRATE/minfreq);  % spacing in points
    GSTATE.FIELDX = fft(GSTATE.FIELDX);
    GSTATE.FIELDX = fastshift(GSTATE.FIELDX,ndfn);
    GSTATE.FIELDX = ifft(GSTATE.FIELDX .* Hf);
    if isy
        GSTATE.FIELDY = fft(GSTATE.FIELDY);
	GSTATE.FIELDY = fastshift(GSTATE.FIELDY,ndfn);
        GSTATE.FIELDY = ifft(GSTATE.FIELDY .* Hf);
    end

end
delfil = evaldelay(ftype,bw*0.5);
GSTATE.DELAY(1,ich) = GSTATE.DELAY(1,ich)+delfil;
if isy
    GSTATE.DELAY(2,ich) = GSTATE.DELAY(2,ich)+delfil;
end

if GSTATE.PRINT

    %%%% PRINT summary

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===            optfilter             ===\n');
    fprintf(fid,'========================================\n\n');
    fprintf(fid,'Optical filtering of ch. #%d\n\n',ich);
    fprintf(fid,'Filter type: %s\n',ftype);
    fprintf(fid,'3-dB bandwidth: %.4f\n',bw);
    fprintf(fid,'Group delay: %.4f bit times\n',delfil);
    if (nfc ~= GSTATE.NCH)
	fprintf(fid,'\nNote: the filter operated on the unique field\n');
    else
	fprintf(fid,'Average energy before filtering: %.4f\n',Eb_before/tpower);
	fprintf(fid,'Average energy after filtering:  %.4f\n',Eb_after/tpower);
	fprintf(fid,'\nNote: other channels left unchanged\n');
    end
    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT
