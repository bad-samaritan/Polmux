function plotfield(pol,ich,flag,col,fil,bw,ord)

%PLOTFIELD Plot the optical field.
%   PLOTFIELD(POL,ICH,FLAG) plots the power and/or phase of channel ICH
%   contained into the global variable GSTATE.FIELDX if POL='x',
%   GSTATE.FIELDY if POL='y', or both if POL='xy'. POL can also be 'tot'
%   and indicates to plot the overall power in the time domain, or the 
%   overall power of the spectrum. Such option works only with FLAG='p---', 
%   'n---', '--p-' or '--n-', see later.
%
%   FLAG is a 4-char string indicating the type of plot:
%
%       FLAG(1)= 'p': power  in the time domain (normalized to the symbol 
%                time). Use 'n' if you want the power normalized to the
%                transmitted peak power of channel ICH.
%       FLAG(2)= 'a': angle (phase) in the same time domain
%       FLAG(3)= 'p': 10*log10(abs(.)^2) of the spectrum (FFT) in the 
%                frequency domain (normalized to the symbol rate). 
%                Use 'n' for a spectrum normalized to the transmitted peak 
%                power of channel ICH. 
%       FLAG(4)= 'a': angle (phase) of the spectrum in the same frequency 
%                domain.
%
%   A specific char of FLAG should be set to '-' to avoid its plot.
%   E.g. Plot only the power in time and the phase of the spectrum:
%           FLAG = 'p--a'
%   E.g. Plot both power and phase in time and frequency:
%           FLAG = 'papa'
%
%   PLOTFIELD(POL,ICH,FLAG,COL,FIL,BW,ORD)
%
%   With this call the specific channel ICH can be temporary extracted   
%   using the filter FIL having 3-dB bandwidth BW and order ORD (see  
%   MYFILTER). Such a requirement is necessary when the fields are combined 
%   into a unique field (see CREATE_FIELD). Otherwise, if ICH=1 and FIL and  
%   BW are not indicated, all the complete unique field is plotted. The
%   option POL='tot' is unactive with FIL present.
%   
%   COL is the color of the lines (see PLOT. E.g. COL='b-','r-x',etc)
%
%   See also PRINTFIELD
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

% The following are my personal settings for figure Position under
% Matlab [width height] (normalized). Change as you like...
gcf12 = [1.6 0.8];    % [1,2] plots
gcf21 = [0.9 1.6];    % [2,1] plots
gcf22 = [1.5 1.5];    % [2,2] plots
gcfhh = [1 1];        % [>2,>2] plots (Units normalized: full screen)

if ~ischar(flag)
    error('FLAG must be a string');
end
if exist('col') ~= 1
    col = 'b';
elseif ~ischar(col)
    error('COL must be a string (E.g. ''b'',''r-o'')');
end
[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 

if (nargin < 7)
    ord = 0;
end

%Lpol = length(pol);
Ltime = 0;
Lfreq = 0;
if (strcmp(flag(1),'p') || strcmp(flag(1),'n')), Ltime=Ltime+1;end;
if strcmp(flag(2),'a'), Ltime=Ltime+1;end;
if (strcmp(flag(3),'p') || strcmp(flag(3),'n')) ...
 || ((flag(3)>'0') && (flag(3) <='9')), Lfreq=Lfreq+1;end;
if strcmp(flag(4),'a'), Lfreq=Lfreq+1;end;

if strcmp(pol,'tot')
    if ~(strcmp(flag,'p---') || strcmp(flag,'n---') || ...
            strcmp(flag,'--p-') || strcmp(flag,'--n-'))
        error('Option ''tot'' works only with flag ''p---'' or ''n---''');
    end
end
Lnot = length(findstr(flag,'-'));
if Ltime+Lfreq+Lnot ~= 4
    error('Wrong FLAG (E.g. FLAG=''p---'',''p-a-'',etc');
end

nrow = max(Ltime,Lfreq);
if (Ltime == 0) || (Lfreq == 0)
    ncol = 1;
else
    ncol = 2;
end
if Ltime > 0
    time = 0:1/GSTATE.NT:GSTATE.NSYMB-1/GSTATE.NT;
end

isy = ~isempty(GSTATE.FIELDY);  
colt{1} = col;
switch pol
    case 'y'
        nend = double(isy);
    case {'x','tot'}
        nend = 1;
    case 'xy'
        nend = 1+isy;
        fscol1 = findstr(col,'-');
        if isempty(fscol1)
            colt{2} = [col,'--'];
        else
            Lcol = length(col);
            if length(fscol1) == 1
                colt{2}(1:fscol1) = col(1:fscol1);
                colt{2}(fscol1+1) = '-';
                colt{2}(fscol1+2:Lcol+1)= col(fscol1+1:Lcol);
            else
                colt{2}(1:fscol1) = col(1:fscol1);
                colt{2}(fscol1+1:Lcol-1) = col(fscol1+2:Lcol);
            end
        end
    otherwise
        error('Wrong string for POL');
end
if exist('OCTAVE_VERSION') == 0
    fsz = get(0,'DefaultFigurePosition');
    if (nrow == 1) && (ncol == 2)
        set(gcf,'Position',[fsz(1)-(gcf12(1)-1)/2*fsz(3),...
                fsz(2)+fsz(4)*(1-gcf12(2)),gcf12(1)*fsz(3),gcf12(2)*fsz(4)]);
    elseif (nrow == 2) && (ncol == 1)
        set(gcf,'Position',[fsz(1)-(gcf21(1)-1)/2*fsz(3),...
                fsz(2)+fsz(4)*(1-gcf21(2)),gcf21(1)*fsz(3),gcf21(2)*fsz(4)]);
    elseif (nrow == 2) && (ncol == 2)
        set(gcf,'Position',[fsz(1)-(gcf22(1)-1)/2*fsz(3),...
                fsz(2)+fsz(4)*(1-gcf22(2)),gcf22(1)*fsz(3),gcf22(2)*fsz(4)]);       
    elseif (nrow > 2) || (ncol > 2)
        set(gcf,'Units','Normalized','Position',[0 0 gcfhh(1) gcfghh(2)]);
    end
end
        
%%%%% READY TO PLOT

Nfft = GSTATE.NSYMB*GSTATE.NT;
for k=1:nend    % cycle over the polarizations
    kcycle = 1;
    if (nfc == GSTATE.NCH) || ((ich == 1) && (exist('fil') ~= 1))
        if strcmp(pol(k),'x')
            fieldt = GSTATE.FIELDX(:,ich);
        elseif strcmp(pol(k),'y')
            fieldt = GSTATE.FIELDY(:,ich);
        else
            fieldt = sqrt(abs(GSTATE.FIELDX(:,ich)).^2+...
                abs(GSTATE.FIELDY(:,ich)).^2);      
        end
        locpow = GSTATE.POWER(ich);                    
        if exist('fil') == 1
            fieldt = ifft(fft(fieldt).*myfilter(fil,GSTATE.FN,bw*0.5,ord));
        end
    else    % temporary extract the channel
        minfreq = GSTATE.FN(2)-GSTATE.FN(1);    
        maxl=max(GSTATE.LAMBDA);
        minl=min(GSTATE.LAMBDA);
        lamc = 2*maxl*minl/(maxl+minl); %central wavelength    
        deltafn = CLIGHT*(1/lamc-1./GSTATE.LAMBDA(ich)); 
        ndfn = round(deltafn./GSTATE.SYMBOLRATE/minfreq);  % spacing in points

        if strcmp(pol(k),'x')
            fieldt = fft(GSTATE.FIELDX);
            fieldt = fastshift(fieldt,ndfn);	% right circular shift
            fieldt = ifft(fieldt.*myfilter(fil,GSTATE.FN,bw*0.5,ord));
        elseif strcmp(pol(k),'y')
            fieldt = fft(GSTATE.FIELDY);      
            fieldt = fastshift(fieldt,ndfn);	% right circular shift
            fieldt = ifft(fieldt.*myfilter(fil,GSTATE.FN,bw*0.5,ord));
        else
            fieldt = fft(GSTATE.FIELDX);
            fieldt = fastshift(fieldt,ndfn);	% right circular shift
            fieldt = ifft(fieldt.*myfilter(fil,GSTATE.FN,bw*0.5,ord));
            fieldty = fft(GSTATE.FIELDY);
            fieldty = fastshift(fieldty,ndfn);	% right circular shift
            fieldty = ifft(fieldty.*myfilter(fil,GSTATE.FN,bw*0.5,ord));
            fieldt = sqrt(abs(fieldt).^2+abs(fieldty).^2);
        end
        locpow = GSTATE.POWER(ich);        
    end
    if Lfreq > 0
        %fieldf = fft(fieldt)/sqrt(Nfft);
        if strcmp(pol,'tot')
            fieldf = sqrt(abs(fft(GSTATE.FIELDY(:,ich))/Nfft).^2 + ...
                abs(fft(GSTATE.FIELDY(:,ich))/Nfft).^2);
        else
            fieldf = fft(fieldt)/Nfft;
        end
    end
   
    if strcmp(flag(1),'p')  % power in the time domain
        if ~(exist('OCTAVE_VERSION','var') || (sscanf(version,'%d')<7))
            ax_v_time(kcycle)=subplot(nrow,ncol,kcycle);
        else
            subplot(nrow,ncol,kcycle);
        end
        hold on;
        plot(time,abs(fieldt).^2,colt{k});
        xlabel('time [symbols]')
        ylabel('power [mW]')
        grid on
        if k > 1
            legend('x','y');
        end
        kcycle = kcycle+1;
        drawnow;
    elseif strcmp(flag(1),'n')  % normalized power in the time domain
        if ~(exist('OCTAVE_VERSION') || (sscanf(version,'%d')<7))
            ax_v_time(kcycle)=subplot(nrow,ncol,kcycle);
        else
            subplot(nrow,ncol,kcycle);
        end
        hold on;
        plot(time,abs(fieldt).^2/locpow,colt{k});
        xlabel('time [symbols]')
        ylabel('normalized power [a.u.]')
        grid on
        if k > 1
            legend('x','y');
        end
        kcycle = kcycle+1;
        drawnow;
    end    
    if strcmp(flag(3),'p')  % power in the frequency domain
        if ~(exist('OCTAVE_VERSION') || (sscanf(version,'%d')<7))
            ax_v_freq(kcycle)=subplot(nrow,ncol,kcycle);
        else
            subplot(nrow,ncol,kcycle);
        end
        hold on;
        plot(fftshift(GSTATE.FN),fftshift(10*log10(abs(fieldf).^2)),colt{k});
        xlabel('frequency [a.u.]')
        ylabel('PSD   [dB]')   
        grid on
        if k > 1
            legend('x','y');
        end        
        kcycle = kcycle+1;
        drawnow;
    elseif strcmp(flag(3),'n')  % normalized power in the frequency domain
        if ~(exist('OCTAVE_VERSION') || (sscanf(version,'%d')<7))        
            ax_v_freq(kcycle)=subplot(nrow,ncol,kcycle);
        else
            subplot(nrow,ncol,kcycle);
        end
        hold on;
        plot(fftshift(GSTATE.FN),fftshift(10*log10(abs(fieldf).^2/locpow)),colt{k});
        xlabel('frequency [a.u.]')
        ylabel('normalized PSD   [dB]')   
        grid on
        if k > 1
            legend('x','y');
        end        
        kcycle = kcycle+1;
        drawnow;
    end    
    if strcmp(flag(2),'a')  % phase in the time domain
        if ~(exist('OCTAVE_VERSION') || (sscanf(version,'%d')<7))
            ax_v_time(kcycle)=subplot(nrow,ncol,kcycle);
        else
            subplot(nrow,ncol,kcycle);       
        end
        hold on;
        plot(time,angle(fieldt),colt{k});
        angles={'-1';'-3/4';'-1/2';'-1/4';'0';'1/4';'1/2';'3/4';'1'};
        set(gca,'YTick',[-pi:pi/4:pi]);
        set(gca,'YTickLabel',angles)
        xlabel('time [symbols]')
        ylabel('phase   [rad/\pi]')    
        grid on
        if k > 1
            legend('x','y');
        end        
        kcycle = kcycle+1;
        drawnow;
    end
    if strcmp(flag(4),'a')  % phase in the frequency domain
        if ~(exist('OCTAVE_VERSION') || (sscanf(version,'%d')<7))
            ax_v_freq(kcycle)=subplot(nrow,ncol,kcycle);
        else
            subplot(nrow,ncol,kcycle);
        end
        hold on;
        plot(fftshift(GSTATE.FN),fftshift(angle(fieldf)),colt{k});
        angles={'-3/4';'-1/2';'-1/4';'0';'1/4';'1/2';'3/4'};
        set(gca,'YTick',[-3*pi/4:pi/4:3*pi/4]);
        set(gca,'YTickLabel',angles)
        xlabel('frequency [a.u.]')
        ylabel('Fourier phase [rad/\pi]')    
        grid on
        if k > 1
            legend('x','y');
        end        
        kcycle = kcycle+1;
        drawnow;
    end
end

if ~(exist('OCTAVE_VERSION') || (sscanf(version,'%d')<7))
    if exist('ax_v_time','var') && length(ax_v_time)>1
        linkaxes(ax_v_time(ax_v_time~=0),'x');
    end
    if exist('ax_v_freq','var') && length(ax_v_freq)>1
        linkaxes(ax_v_freq(ax_v_freq~=0),'x');
    end
end

