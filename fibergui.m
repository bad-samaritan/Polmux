function fibergui(Lf,alphadB,Aeff,n2,lambda,Dc,Slope,dzmax,...
    dphimax,flag,infoax)

%FIBERGUI Optical fiber in the nonlinear regime (GUI tool).
%   FIBERGUI(LF,ALPHAdB,AEFF,N2,LAMBDA,DC,SLOPE,DZMAX,DPHIMAX,FLAG,INFOAX) 
%   solves the nonlinear Schroedinger equation (NLSE) in absence of 
%   polarization  effects using a GUI interface. The fiber parameters are:
%
%       LF:         length [m]
%       ALPHAdB:    attenuation [dB/km]
%       Aeff:       effective area [um^2]
%       n2:         nonlinear index [m^2/W]
%       LAMBDA:     wavelength [nm] of DC 
%       DC:         fiber chromatic dispersion coefficient [ps/nm/km]
% 					@ LAMBDA
%       SLOPE:      fiber slope, i.e. derivative of DC [ps/nm^2/km]
%       DZMAX:      max. step for the split-step algorithm [m]
%       DPHIMAX:    max. nonlinear phase rotation in each step [rad]
%
%   The NLSE is solved by a split-step Fourier algorithm with a variable
%   step so as to have a maximum nonlinear phase rotation into each step
%   equal to DPHIMAX. However, the step cannot be larger than DZMAX.
%
%   FLAG is a string of three characters governing the type of propagation. 
%   The first character is 'g' if GVD (i.e. beta2,beta3) is on or '-' in 
%   absence of GVD. Note that with 'sepfields' in CREATE_FIELD this 
%	function accounts for the walkoff effect even with the GVD flag set to 
%	'-'. The third is 's' if SPM is on or '-' in absence of SPM. 
%   Likewise, the fourth character is 'x' or '-' in presence/absence of 
%   XPM. The most complete case is FLAG='gsx' and corresponds to 
%   propagation in presence of fiber GVD+SPM+XPM. 
%
%   E.g. Propagation with GVD+XPM, without SPM -> FLAG='g-x'
%
%   The third character of FLAG is active only with channels separated 
%   (see option 'sepfields' in CREATE_FIELD). In this case, the propagation 
%   neglects the effect of four-wave mixing, which can be taken in account 
%   only by combining all channels into a unique field and hence it is a 
%   special case of SPM.
%
%   INFOAX governs the GUI tool. The options of INFOAX are:
%
%       INFOAX.ch = channel under investigation;
%       INFOAX.flag1d = [a b c]. a=1 -> plot the power. b=1 -> plot the
%           phase. c=1 -> plot the chirp. a or b or c=0 -> don't plot.
%       INFOAX.flag3d = same form as INFOAX.flag1d, but for 3D plot
%       INFOAX.ch = channel to be plotted. E.g. INFOAX.ch = [1 3] indicates 
%           that only channels 1 and 3 will be plotted.
%       INFOAX.axprop = cell array containing valid pairs of axes
%           properties (see AXES for more details).
%           E.g. INFOAX.axprop = {'XLim',[0 10],'YLim',[0 2],'FontSize',24}
%
%   Note 1: FIBERGUI does not account for PMD yet. In future versions of 
%		Optilux it will be implemented like FIBER.
%	Note 2: This function contains some known minor bugs.
%	Note 3: FIBERGUI does not work under Octave.
%
%   See Also FIBER, FIELDGUI, UPDSEL
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
MAXSTEP = 131072;
DEF_MAX3D = 30;

[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
Nfft = GSTATE.NSYMB*GSTATE.NT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK AND INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 11
    error('missing some input arguments');
end 
fls = [0 0 0];  % flag converted into array of numbers
switch lower(flag)
    case 'g--'  % only GVD
        fls(1) = 1;
        dphimaxt = Inf;
        dzmaxt = Lf;
    case '-s-'  % only SPM
        fls(2) = 1;
        dphimaxt = Inf;
        dzmaxt = Lf;      
    case '--x'  % only XPM 
        if nfc == 1
            error(['flag ''--x'' available only for',...
                'channels separated']);
        end
        fls(3) = 1;
        dphimaxt = Inf;
        dzmaxt = Lf;   
    case 'gs-'  % GVD+SPM
        fls(1:2) = 1;
        dphimaxt = dphimax;
        dzmaxt = dzmax;              
    case 'g-x'  % GVD+XPM
        if nfc == 1
            error(['flag ''g-x'' available only for',...
                'channels separated']);
        end        
        fls([1,1],[1,3]) = 1;
        dphimaxt = dphimax;
        dzmaxt = dzmax;              
    case '-sx'  % SPM+XPM
        fls(2) = 1;
        if (nfc ~= 1), fls(3) = 1;end;
        dphimaxt = Inf;
        dzmaxt = Lf;      
    case 'gsx'  % GVD+SPM+XPM: complete case
        fls(1:2) = 1;
        if (nfc ~= 1), fls(3) = 1;end;        
        dphimaxt = dphimax;
        dzmaxt = dzmax;              
    otherwise
        error(['wrong flag. E.g. ',...
            'flag can be ''g--'',''gs-'',''-s-'',etc']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alphalin = (log(10)*1e-4)*alphadB;      % [m^-1]
if alphalin == 0
    Leff = Lf;      % effective length [m]
else
    Leff = (1-exp(-alphalin*Lf))/alphalin;
end
b20 = -lambda^2/2/pi/CLIGHT*Dc*1e-6; % beta2 [ns^2/m] @ lambda
b30 = (lambda/2/pi/CLIGHT)^2*(2*lambda*Dc+lambda^2*Slope)*1e-6;  
                                     % beta3 [ns^3/m] @ lambda

maxl = max(GSTATE.LAMBDA);
minl = min(GSTATE.LAMBDA); 
lamc = 2*maxl*minl/(maxl+minl); %central wavelength

% Domega_ik: [1/ns]. "i" -> at ch. i, "0" -> at lambda, "c" -> at lamc
Domega_i0 = 2*pi*CLIGHT*(1./GSTATE.LAMBDA-1/lambda);    
Domega_ic = 2*pi*CLIGHT*(1./GSTATE.LAMBDA-1/lamc);    
Domega_c0 = 2*pi*CLIGHT*(1./lamc-1/lambda); 
b1 = b20*Domega_ic+0.5*b30*(Domega_i0.^2-Domega_c0.^2);  %ch's beta1 [ns/m]
if nfc == 1
    beta1 = 0;    % [ns/m] @ lamc
    Domega_i0 = 2*pi*CLIGHT*(1./lamc-1/lambda);        
else                                     
    beta1 = b1;   % [ns/m] @ GSTATE.LAMBDA
end
beta2 = b20+b30*Domega_i0;  % beta2 [ns^2/m]@ lamc (nfc=1) or GSTATE.LAMBDA
                            % (nfc ~= 1).
Dch = Dc+Slope*(GSTATE.LAMBDA-lambda);  % dispersion of the channels
% Note: channel 1 is the reference channel. All the other channels shift
%   (beta1 effect) with respect to it.

Ld = Inf*ones(1,GSTATE.NCH);
fzb2 = find(beta2);
Ld(fzb2) = 1./(GSTATE.SYMBOLRATE^2*abs(beta2(fzb2))); % dispersion length [m]
if b30 ~= 0
    Lds = 1./(GSTATE.SYMBOLRATE^3*abs(b30));      % slope length [m]
else
    Lds = Inf;
end

betat = zeros(Nfft,nfc);    % beta, whose Taylor coeffs are the betai
omega = 2*pi*GSTATE.SYMBOLRATE*GSTATE.FN';     % angular frequency [rad/ns]

for kch=1:nfc
    betat(:,kch) = omega*beta1(kch)+0.5*omega.^2*beta2(kch)+...
        omega.^3*b30/6;   % beta coefficient [1/m]
end
% gamma nonlinear index [1/mW/m]
gam = 2*pi*n2./(GSTATE.LAMBDA*Aeff)*1e18;
Lnl = 1./(gam.*GSTATE.POWER);    % nonlinear length [m]
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% UPDATE DELAY AND DISP %%%%%%%%%%%%%%%%%%%%%%%%%%%

loc_delay = Lf*GSTATE.SYMBOLRATE.*b1;
GSTATE.DELAY = GSTATE.DELAY+loc_delay;
GSTATE.DISP = GSTATE.DISP + Dch*Lf*1e-3;    % cumulated dispersion [ps/nm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SSFM ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(infoax,'comp')
    infoax.comp = 0;
end
if (sum(infoax.flag3d)) 
    infoax.y3D = ones(1,Nfft);
    if ~isfield(infoax,'max3d'), infoax.max3d = DEF_MAX3D; end;
end
selector = 2;   % default gui operation: play
time = 0:1/GSTATE.NT:GSTATE.NSYMB-1/GSTATE.NT;   % x-axis for the plot
infoax=fieldgui(infoax);  %  initialize gui interactive plot
allstep = zeros(1,MAXSTEP); % vector of all tested steps [m]
ncycle = 1;      % number of cycles
dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,GSTATE.FIELDX,nfc);
halfalpha = 0.5*alphalin; 
allstep(1) = dz;
firstdz = dz;
zprop = dz;             % running distance [m]
while zprop < Lf        % all steps except the last
    GSTATE.FIELDX = nl_step(alphalin,gam,dz,GSTATE.FIELDX,...
        nfc,fls(2),fls(3));       % 1/3) NON-linear step
    
    if fls(1) == 1  %    %%% GVD %%%
        GSTATE.FIELDX = lin_step(betat,dz,GSTATE.FIELDX); % 2/3) linear step    
    end
    GSTATE.FIELDX = GSTATE.FIELDX*exp(-halfalpha*dz);    % 3/3) attenuation
    %%%%% GUI
    infoax.dist = zprop;    % running distance for the title of the figure    
    if infoax.comp == 1
        opt = exp(halfalpha*zprop)*...
            lin_step(-betat,zprop,GSTATE.FIELDX(:,infoax.ch));
    else
        opt  = GSTATE.FIELDX(:,infoax.ch)*exp(0.5*alphalin*zprop); % gui field        
    end
    actplot(time,opt,infoax,sign(dz));
    if selector ~= 2    % if not play
        uiwait; % wait for a button pressure
    end
    selector = get(gcf,'UserData');    % get button values
    switch selector  % operation to do
        
        %%%%%%
        case 1       % REWIND
        %%%%%%    
            if (ncycle >= 1) 
                dz = -allstep(ncycle);  % go backward
                ncycle = ncycle - 1;
            else
                dz = 0;
            end
            
        %%%%%%    
        case 2       % PLAY
        %%%%%%    
            twait = get(infoax.slider,'Value');  % get the pause
             dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,GSTATE.FIELDX,...
                nfc);  %  go ahead
            pause(twait);
            ncycle = ncycle + 1;
            allstep(ncycle) = dz;            
            
        %%%%%%    
        case 3       % FORWARD
        %%%%%%
            dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,GSTATE.FIELDX,...
                nfc);  %  go ahead
            ncycle = ncycle+1;
            allstep(ncycle)=  dz;            
        
        %%%%%%
        case 4       % PAUSE
        %%%%%%
            dz = 0;            
    end
            
    zprop = zprop+dz;
end 
last_step = Lf-zprop+dz;
infoax.dist = Lf;   
GSTATE.FIELDX = nl_step(alphalin,gam,last_step,GSTATE.FIELDX,...
    nfc,fls(2),fls(3)); % NON-linear step
if fls(1) == 1  %    %%% GVD %%%
    GSTATE.FIELDX = lin_step(betat,last_step,GSTATE.FIELDX); % linear step    
end
GSTATE.FIELDX = GSTATE.FIELDX*exp(-halfalpha*last_step);    % attenuation
opt  = GSTATE.FIELDX(:,infoax.ch)*exp(0.5*alphalin*zprop);
actplot(time,opt,infoax,sign(dz));

if GSTATE.PRINT
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% PRINT SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===       Interactive fiber          ===\n');
    fprintf(fid,'========================================\n\n');
    fprintf(fid,'Fiber parameters:\n\n');
    fprintf(fid,'Length:%17.3f  [km]\n',Lf*1e-3);
    fprintf(fid,'Attenuation:%12.2f  [dB/km] (Leff = %7.3f [km])\n',...
	alphadB,Leff*1e-3);
    fprintf(fid,'lambda of Dc:%11.2f  [nm]\n',lambda);
    fprintf(fid,'Dc:%21.4f  [ps/nm/km]\n',Dc);
    fprintf(fid,'Slope:%18.4f  [ps/nm^2/km]\n',Slope);
    fprintf(fid,'n2:%21.2e  [m^2/W]\n',n2);
    fprintf(fid,'Aeff:%19.2f  [um^2]\n\n',Aeff);
    fprintf(fid,'Propagation type: ''%s''\n\n',flag);
    fprintf(fid,'Max NL phase rotation x step: %-6.2g  [rad]\n',dphimax);
    fprintf(fid,'Max step: %.2e  [m]\n',dzmax);
    fprintf(fid,'Initial step: %.2e  (num. steps: %d)\n\n',firstdz,ncycle);
    fprintf(fid,'Channel properties (Ld: disp. length. Lnl: NL length):\n\n');
    for kch=1:GSTATE.NCH
	fprintf(fid,'ch. #%.2d: Dc = %.4f  [ps/nm/km]  (Ld = %3.2e [km])\n',...
            kch,Dch(kch),Ld(kch)*1e-3);
	fprintf(fid,'\t gamma = %.3e [1/mW/km] (Lnl = %3.2e [km])\n',...
            gam(kch)*1e3,Lnl(kch)*1e-3);
	fprintf(fid,'\t sqrt(Ld/Lnl) = %.4f\n',sqrt(Ld(kch)/Lnl(kch)));
    end
    fprintf(fid,'\nSlope length Lds: %-3.2e  [km]\n',Lds*1e-3);

    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT

%-----------------------------------------------------
function dz=nextstep(dzmax,phimax,gam,alphalin,u,nfc)

%NEXTSTEP step for the SSFM algorithm
%   DZ=NEXTSTEP(DZMAX,PHIMAX,GAM,ALPHALIN,U,NFC) evaluates the step of the 
%   SSFM algorithm for a fiber having nonlinear gamma coefficient GAM 
%   [1/mW/km]and attenuation ALPHALIN [m^-1]. U is the matrix of electrical
%   fields. NFC is the number of channels (equal to 1 if all channels are 
%   combined into a unique field). The step corresponds to a maximum 
%   nonlinear phase rotation equal to PHIMAX, under the constraint that it 
%   cannot be greater than DZMAX.

Umax = max(real(u).^2+imag(u).^2);               % largest normalized power
Pmax = max(gam(1:nfc).*Umax);   % largest gamma*power
leff = phimax/Pmax;            % effective length of the step
dl = alphalin*leff;                  % ratio effective length/attenuation length

if dl >= 1
    dz = dzmax;
else
    if alphalin == 0
        step = leff;
    else
        step = -1/alphalin*log(1-dl);
    end
    if step > dzmax
        dz = dzmax;
    else
        dz = step;
    end
end

%-----------------------------------------------------------
function u = lin_step(beta,dz,u)

%LIN_STEP linear fiber
%   U=LIN_STEP(BETA,DZ,U) propagates the electric field U
%   into a purely linear fiber having beta factor BETA, length DZ. 


Hf = fastexp(-beta*dz); %  Fast exponential: fastexp(x) = exp(i*x)

u = ifft( fft(u) .* Hf);

%-----------------------------------------------------------
function u = nl_step(alphalin,gam,dz,u,nfc,spm,xpm)

%NL_STEP nonlinear fiber
%   U=NL_STEP(ALPHALIN,GAM,DZ,U,NFC,SPM,XPM) propagates the electric field 
%   U into a purely nonlinear fiber having a nonlinear fiber having
%   nonlinear coefficient GAM [1/mW/km], attenuation ALPHA_LIN [m^-1]. NFC 
%   is the number of channels. SPM and XPM are 1 or 0 if self-phase 
%   modulation and cross-phase modulation have to be taken in account, 
%   respectively.


if alphalin == 0
    leff = dz;
else
    leff = (1-exp(-alphalin*dz))/alphalin;  % effective length of dz
end
for kch=1:nfc
    if spm          %  Fast exponential: fastexp(x) = exp(i*x)
        u(:,kch) = u(:,kch) .* fastexp(-gam(kch)*(real(u(:,kch)).^2+...
            imag(u(:,kch)).^2)*leff);
    end
    if xpm
        n = [1:kch-1,kch+1:nfc];    % XPM channels
        for nch=1:nfc-1
            u(:,kch) = u(:,kch) .* ...
                fastexp(-2*gam(n(nch))*(real(u(:,n(nch))).^2+...
                imag(u(:,n(nch))).^2)*leff);
        end
    end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function actplot(x,y,infoax,yout)

%ACTPLOT update plot
%   ACTPLOT(X,Y,INFOAX,YOUT) updates the plot Y=f(X) or Y=f(X,k), being k a
%   counter for 3D plot. INFOAX contains the plot properties. YOUT is the
%   type of step (+/-1).
%
%   Author: Paolo Serena, 2007
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2007  Paolo Serena, <serena@tlc.unipr.it>
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


set(infoax.fig,'Name',['Interactive Fiber. Distance = ',...
    num2str(infoax.dist),' m']);
n = 1;
ksig = 1;
while n <= infoax.numaxes
    if infoax.flag1d(1)    %%% Power
        plot1dim(x,abs(y(:,ksig)).^2,infoax.axes(n));
        n = n+1;
    end
    if infoax.flag1d(2)    %%% Phase
        plot1dim(x,angle(y(:,ksig)),infoax.axes(n));
        n = n+1;
    end
    if infoax.flag1d(3)    %%% Chir
        plot1dim(x,angle(y(:,ksig)),infoax.axes(n));        
        n = n+1;
    end
    if infoax.flag3d(1)    %%% Power (3D)
        plot3dim(x,infoax.y3D,abs(y(:,ksig)).^2,infoax.axes(n),yout,...
            infoax.max3d);
        n = n+1;
    end
    if infoax.flag3d(2)    %%% Phase (3D)
        plot3dim(x,infoax.y3D,angle(y(:,ksig)),infoax.axes(n),yout,...
            infoax.max3d);        
        n = n+1;
    end
    if infoax.flag3d(3)    %%% Chirp (3D)
        plot3dim(x,infoax.y3D,angle(y(:,ksig)),infoax.axes(n),yout,...
            infoax.max3d);                
        n = n+1;
    end    
    ksig = ksig+1;
end


%------------------------------------------------------------
function plot1dim(x,y,axn)

%PLOT1DIM(X,Y,AXN) 1D interactive plot.
%   X: x-axis. Y:y-axis. AXN: axes. Function of fieldgui.m

axes(axn);
oldline = findobj(gca,'Type','line');
delete(oldline);
line(x,y,'EraseMode','normal');


%------------------------------------------------------------
function plot3dim(x,y,z,axn,yout,max3d)

%PLOT3DIM(X,Y,AXN) 3D interactive plot.
%   X: x-axis. Y:y-axis. Z.z-axis. AXN: axes. YOUT= type of step.
%   MAX3D: max number of lines before clearing.
%   Function of fieldgui.m

axes(axn);
oldline = findobj(gca,'Type','line');
LO = length(oldline);
if (yout <= 0) % backward
    if (LO) delete(oldline(1)); end;
else           % play or forward 
    if (LO >= max3d) 
        delete(oldline); 
        LO = 0;
    end;
    line(x,(LO+1)*y,z,'EraseMode','normal');
end
