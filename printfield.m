function printfield(pol,ch,name,flag,ncycle,fil,bw,ord)

%PRINTFIELD Print the optical field to file.
%   PRINTFIELD(POL,ICH,NAME,FLAG) prints the electric field of channel ICH
%   into the file NAME in the global directory GSTATE.DIR. 
%   POL can be either 'x' for x-polarization, 'y' for y-polarization or 
%   'xy' for printing both polarizations.
%
%   FLAG is a string that indicates the type of files to be printed. It can
%   be 'papa' in the most complete case, indicating that it will be printed
%   both the power and angle (phase) in the time domain (first two chars 
%   of FLAG) and the power (abs(.)^2) and phase of the spectrum (last two 
%   chars of FLAG). If you don't want to print a specific component set it 
%   to '-'.
%
%   E.g. Print only the power in time and the phase of the spectrum:
%           FLAG = 'p--a'
%
%   PRINTFIELD(POL,ICH,NAME,FLAG,NCYCLE,FIL,BW)
%
%   NCYCLE is a counter that will be added to the output file name, and can 
%   be useful in presence of an external cycle to the function.
%
%   If all channels are combined into a unique field, you can temporary
%   extract the ICH channel with the filter FIL having 3-dB bandwidth BW 
%   (see MYFILTER). In this case and in absence of an external cycle, leave
%   NCYCLE=[]. However, if only the spectrum is requested in FLAG, ICH is
%   ignored and all the complete spectrum is saved into one file.
%
%   PRINTFIELD(POL,ICH,NAME,FLAG,NCYCLE,FIL,BW,ORD) use the special 
%   parameter ORD for special filters. E.g. ORD is the supergauss order for
%   the supergaussian filter (see MYFILTER).
%
%   The power will be printed in a directory with the suffix ".MOD",
%   while the phase with the suffix ".ANG".
%
%   E.g. POL='x',CH=1,NCYCLE=81,FLAG='p---',NAME='Tx' 
%               -> creates tempx_ch01_Tx_Cycle_0081
%
%   Note: PRINTFIELD presumes that the user initialized GSTATE.DIR in
%   RESET_ALL. If not, a temporary directory is used and a warning is
%   displayed.
%
%   See also PLOTFIELD, RESET_ALL
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

MAXD = 32768;   % maximum number of points for using fprintf
ncount = 1;     % counter used in the "print summary" section

if length(flag) ~= 4
    error('FLAG must be of length 4');
end
if (length(findstr(flag,'p'))+length(findstr(flag,'a'))+...
        length(findstr(flag,'-'))) ~= 4
    error('Wrong FLAG (E.g. FLAG=''p---'',''p-a-'',etc');
end
Nfft = GSTATE.NSYMB*GSTATE.NT;
[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
                                    % for all channels.
if strcmp(flag(1),'-') && strcmp(flag(2),'-') && (nfc ~= GSTATE.NCH)    
    ich = 1;    % extract the overall spectrum
else
    ich = ch;
end

isy = ~isempty(GSTATE.FIELDY);     
switch pol
    case 'y'
        nend = double(isy);
    case 'x'
        nend = 1;
    case 'xy'
        nend = 1+isy;
    otherwise
        error('Wrong string for POL');
end
if Nfft > MAXD
    types = 1;  % use save for printing (faster than fprintf)
else
    types = 0;  % use fprintf (smaller file and more readable)
end


%%%%%%%%%%%%%%%%%
for n=1:nend       % POLARIZATIONS CYCLE
%%%%%%%%%%%%%%%%%

%%%% EXTRACT THE CHANNEL

if (nfc == GSTATE.NCH) || (ich == 1)
    if strcmp(pol(n),'x')
        field = GSTATE.FIELDX(:,ich);
    else
        field = GSTATE.FIELDY(:,ich);
    end
else    % temporary extract the channel
    minfreq = GSTATE.FN(2)-GSTATE.FN(1); 
    maxl=max(GSTATE.LAMBDA);
    minl=min(GSTATE.LAMBDA);
    lamc = 2*maxl*minl/(maxl+minl); %central wavelength   
    deltafn = CLIGHT*(1/lamc-1./GSTATE.LAMBDA(ich)); % frequency spacing
    ndfn = round(deltafn./GSTATE.SYMBOLRATE/minfreq);  % spacing in points
    if (nargin < 8)
        ord = 0;
    end 
    if strcmp(pol(n),'x')
        field = fft(GSTATE.FIELDX);     
        field = fastshift(field,ndfn);   
        field = ifft(field.*myfilter(fil,GSTATE.FN,bw*0.5,ord));
    else
        field = fft(GSTATE.FIELDY);
        field = fastshift(field,ndfn);
        field = ifft(field.*myfilter(fil,GSTATE.FN,bw*0.5,ord));
    end    
end

%%%% SET DIRECTORIES
if ~isfield(GSTATE,'DIR')
    GSTATE.DIR = tempname('.');
    mkdir(GSTATE.DIR);
    warning('optilux:printfield','Missing output directory in RESET_ALL. %s used',...
        GSTATE.DIR);
end

dirp = [GSTATE.DIR,'/',GSTATE.DIR,'.MOD/'];    
dira = [GSTATE.DIR,'/',GSTATE.DIR,'.ANG/'];    

%%%% SET LABELS
chlab = num2str(ich,'%.2d');    %if ich < 10, adds a '0'

if exist('ncycle') == 1
    if ~isempty(ncycle)
        cyclab = ['_Cycle_',num2str(ncycle,'%.4d')];    
                                % if an external cycle is present, sign it!
    else
        cyclab = '';
    end
else
    cyclab = '';
end

if (flag(1) == 'p') || (flag(2) == 'a')
    st = 1/GSTATE.NT;              % sampling time (norm. to bit time)
    time = 0:st:GSTATE.NSYMB-st;    % time axis
end

%%%% WRITE TO FILE

if flag(1) == 'p'   % TIME-POWER
    nomewrite(ncount,:) = [dirp,'temp',pol(n),'_ch',...
            chlab,'_',name,cyclab,'.dat'];
    printallf(nomewrite(ncount,:),time,abs(field).^2,types);
    ncount = ncount+1;
end
if flag(2) == 'a'   % TIME-ANGLE
    nomewrite(ncount,:) = [dira,'temp',pol(n),'_ch',...
            chlab,'_',name,cyclab,'.dat'];
    printallf(nomewrite(ncount,:),time,angle(field),types);
    ncount = ncount+1;
end
if flag(3) == 'p'   % FREQ-POWER
    nomewrite(ncount,:) = [dirp,'freq',pol(n),'_ch',...
            chlab,'_',name,cyclab,'.dat'];
    printallf(nomewrite(ncount,:),fftshift(GSTATE.FN),...
        fftshift(abs(fft(field)).^2),types);
    ncount = ncount+1;
end
if flag(4) == 'a'   % FREQ-ANGLE
    nomewrite(ncount,:) = [dira,'freq',pol(n),'_ch',...
            chlab,'_',name,cyclab,'.dat'];
    printallf(nomewrite(ncount,:),fftshift(GSTATE.FN),...
        fftshift(angle(fft(field)).^2),types);
    ncount = ncount+1;
end   

%%%
end % end for n=nstart:nend
%%%

if GSTATE.PRINT

    %%%% PRINT summary

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===            printfield            ===\n');
    fprintf(fid,'========================================\n\n');
    for k=1:ncount-1
	fprintf(fid,'%s\n',nomewrite(k,:));
    end
    fprintf(fid,'\nPrinted to file.\n');
    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
    
end  % if GSTATE.PRINT    



