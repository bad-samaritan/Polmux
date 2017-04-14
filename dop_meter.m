function varargout=dop_meter(ich,nfig_flag,fil,bw,ord)

%DOP_METER Computes the Degree Of Polarization of the Optical field.
%   DOP=DOP_METER(ICH) returns the Degree Of Polarization (DOP) of
%   partially polarized light. ICH is the channel number. 
%
%   DOP=DOP_METER(ICH,NFIG_FLAG) with NFIG_FLAG~=0 or nonempty indicates 
%   that the states of polarization (SOP) of channel ICH will be plotted in 
%   the current figure. NFIG_FLAG can be a char, e.g. NFIG_FLAG='b', which 
%   means that the SOP plot will use that color for each sample. Otherwise, 
%   NFIG_FLAG can be a struct with fields:
%
%   NFIG_FLAG.col = color of plot (e.g. 'b'). If equal to the empty string,
%       only the average DOP is plotted.
%   NFIG_FLAG.avgcol = color of the average DOP (e.g. 'r')
%   NFIG_FLAG.norm = 'y': The SOP for each point having power greater than 
%       the signal average power is plotted with modulus 1, and hence lies 
%       on the radius 1 sphere. The other points are discarded.
%   NFIG_FLAG.view = 'best': The best view of the DOP on the sphere.
%
%   The SOP shows the state of partially polarized light as a cloud of 
%   colored points (one for each time sample) surrounding the "polarized 
%   component" of the field, shown with a black vector (default). The 
%   Poincaré sphere is plotted as a reference.
%   As extreme case, the black circle is in the origin for unpolarized
%   light (natural, light, ASE noise...), since the polarized component is
%   null, while the colored dots all collapse on the black circle for
%   totally polarized light, since there is no unpolarized component.
%
%   DOP=DOP_METER(ICH,NFIG_FLAG,FIL,BW,ORD) works for a unique field (see
%   CREATE_FIELD) and temporary extracts channel ICH with an optical filter
%   FIL of bandwidth BW and optional order ORD (see MYFILTER). 
%
%   [DOP,PHI]=DOP_METER(ICH,NFIG_FLAG,FIL,BW) also returns in PHI.azi and
%   PHI.ell the azimuth and ellipticity of the average SOP.
%
%   See also: POL_SCRAMBLER, SET_SOP, POLARIZER
%
%   Author: Armando Vannucci, 2009
%   Contributed by Paolo Serena, 2009
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
global CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fieldx = GSTATE.FIELDX; % copy
fieldy = GSTATE.FIELDY;

[nfr,nfc] = size(fieldx);    % if nfc == 1 there is only one field 
if exist('nfig_flag','var')
    docol = true;
    if ischar(nfig_flag)
        okplot = 1;
        mycol = [nfig_flag(isletter(nfig_flag)),'.'];
        mycolavg = 'k';
    elseif isstruct(nfig_flag)
        okplot = 1;
        if isfield(nfig_flag,'col')
            coltmp = nfig_flag.col(isletter(nfig_flag.col));
            if strcmp(coltmp,'')
                docol = false;
            else
                mycol = [coltmp,'.'];
            end
        else
            mycol = 'm.';
        end
        if isfield(nfig_flag,'avgcol')
            mycolavg = nfig_flag.avgcol(isletter(nfig_flag.avgcol));
        else
            mycolavg = 'k';
        end
    elseif isempty(nfig_flag) || (nfig_flag == 0), 
        okplot = 0;
    else
        okplot = 1;
        mycol = 'm.';
        mycolavg = 'k';
    end      
else
    okplot = 0;
end
if isempty(fieldy),
    warning('optilux:dop_meter','There is no Y field component; unable to compute DOP');
end
maxl=max(GSTATE.LAMBDA);
minl=min(GSTATE.LAMBDA);
lamc = 2*maxl*minl/(maxl+minl); %central wavelength: 1/lamc = 0.5(1/maxl+1/minl)
if nfc ~= GSTATE.NCH % temporary extract channel ich
    if (exist('fil','var') == 0) || (exist('bw','var') == 0)
        error('missing some optical filter parameters for temp extraction');
    end
    if ~exist('ord','var')
        ord = 0;     % not using special filter
    end
    minfreq = GSTATE.FN(2)-GSTATE.FN(1);
    deltafn = CONSTANTS.CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich)); % frequency spacing
    ndfn = round(deltafn./GSTATE.SYMBOLRATE/minfreq);  % spacing in points
    fieldx = fft(fieldx);
    fieldx = fastshift(fieldx,ndfn);	% right circular shift
    fieldx = ifft(fieldx.*myfilter(fil,GSTATE.FN,bw*0.5,ord));
    
    fieldy = fft(fieldy);
    fieldy = fastshift(fieldy,ndfn);	% right circular shift
    fieldy = ifft(fieldy.*myfilter(fil,GSTATE.FN,bw*0.5,ord));
    nch = 1;
else
    nch = ich;
end   
%%%%%%%%%%%%%%%%%%%%%%%% EVALUATION + VISUALIZATION OF DOP %%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVALUATION OF DOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate Stokes coordinates, for every sample and every i-th channel
S0 = abs(fieldx(:,nch)).^2 + abs(fieldy(:,nch)).^2;
S1 = abs(fieldx(:,nch)).^2 - abs(fieldy(:,nch)).^2;
S2 =  2.*real( fieldx(:,nch).*conj(fieldy(:,nch)));
S3 = -2.*imag( fieldx(:,nch).*conj(fieldy(:,nch)));
% Evaluate time-averaged Stokes coordinates and DOP, for every i-th channel
S0avg = mean(S0);
S1avg = mean(S1);
S2avg = mean(S2);
S3avg = mean(S3);
% Evaluation of the Output argument of the function
DOP = sqrt(S1avg^2+S2avg^2+S3avg^2)/S0avg;

varargout{1} = DOP;
if nargout == 2
    phi.azi= sign(S2).*acos(S1./sqrt(S0.^2-S3.^2))/2 + (S2==0).*(1-sign(S1))*pi/4;
    phi.ell= asin(S3./S0)/2;
    phi.aziavg= sign(S2avg)*acos(S1avg/sqrt(S0avg^2-S3avg^2))/2 + ...
        (S2avg==0).*(1-sign(S1avg))*pi/4;
    phi.ellavg= asin(S3avg/S0avg)/2;

    varargout{2} = phi;
end
%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION ON SPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%

if okplot
    if isempty(get(gca,'ZTickLabel'))
        [Xsph,Ysph,Zsph]= sphere;
        mesh(Xsph,Ysph,Zsph,'FaceAlpha',0.6);
        axis equal;
        axis ([-1 1 -1 1 -1 1]);
        %         title(['Polarization states of channel ',num2str(ich),': ',num2str(Nfft),' time samples (magenta)',...
        %             'Polarized component (black) with power P_{pol}. Sphere radius rescaled to mean power']);
        xlabel('s1'); ylabel('s2');zlabel('s3');

        %         lightangle(90,45);
        colormap summer;
        hold on;
    end

    if isfield(nfig_flag,'view') && strcmp(nfig_flag.view,'best')
        view(phi.aziavg*180/pi*2+90,phi.ellavg*180/pi*2);
    else
        view(110,20); % default
    end
    if isfield(nfig_flag,'norm') && strcmp(nfig_flag.norm,'y')
        sphRN = S0;
        iiso = S0 < S0avg;
        S1(iiso) = NaN;  % do not plot these low power points
        S2(iiso) = NaN;
        S3(iiso) = NaN;
    else
        sphRN = S0avg;    % Sphere radius= mean power;
    end

    sphRadius = S0avg;    % Sphere radius= mean power;
    if docol, plot3(S1./sphRN,S2./sphRN,S3./sphRN,mycol);end
    plot3([0, S1avg/sphRadius],[0, S2avg/sphRadius],[0, S3avg/sphRadius],...
        [mycolavg,'-'],'LineWidth',3)
    plot3(S1avg/sphRadius,S2avg/sphRadius,S3avg/sphRadius,[mycolavg,'o'],...
        'LineWidth',5,'MarkerSize',8);
    % NOTE: the lower the DOP,
end;    % if (nfig_flag ~= 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRINT SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GSTATE.PRINT
    outfile = [GSTATE.DIR,'/simul_out'];
    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===             dop_meter             ===\n');
    fprintf(fid,'========================================\n');
    fprintf(fid,'DOP of channel %d \n\n',ich);
    fprintf(fid,'avg. power P=%.3f [mW] =Ppol+Punp; \n',S0avg);
    fprintf(fid,['   polarized component: Ppol=%.3f [mW];\n   ',...
        'unpolarized component: Punp=%.3f [mW]. \n\n'],...
        DOP*S0avg, (1-DOP)*S0avg);
    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT

