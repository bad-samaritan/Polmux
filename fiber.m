function zbrf=fiber(x,flag)

%FIBER Optical fiber in the nonlinear regime.
%   FIBER(X,FLAG) solves the nonlinear Schroedinger equation (NLSE) in 
%   absence of polarization effects, or the Coupled-NLSE (CNLSE) with 
%   polarization effects. X is a structure of fields:
%
%       X.length:     fiber length [m]
%       X.alphadB:    fiber attenuation [dB/km]
%       X.aeff:       fiber effective area [um^2]
%       X.n2:         fiber nonlinear index [m^2/W]
%       X.lambda:     wavelength [nm] of X.disp 
%       X.disp:       fiber chromatic dispersion coefficient [ps/nm/km]
%                     @X.lambda
%       X.slope:      fiber slope, i.e. derivative of X.disp [ps/nm^2/km]
%       X.dzmax:      max. step for the split-step algorithm [m]
%       X.dphimax:    max. nonlinear phase rotation in each step [rad]
%
%   For the solution of the CNLSE, i.e. with two polarizations, there are
%   also the following additional parameters:
%
%       X.dgd:        fiber average differential group delay [symbols]
%       X.nplates:    number of waveplates or trunks for PMD emulation
%       X.manakov:    'yes': Solve the Manakov equation. 'no': Solve the
%                     CNLSE. Default: 'no'.
%
%   In the general case with two polarizations the fiber is the 
%   concatenation of randomly oriented polarization mantaining fibers 
%   (PMF). The user can force the use of a single PMF by adding the 
%   following optional parameters:
%
%       X.db0:        birefringence of the PMF fiber at GSTATE.FN=0 
%       X.theta:      azimuth [rad] of the PMF fiber
%       X.epsilon:    ellipticity [rad] of the PMF fiber
%
%   The NLSE is solved by a split-step Fourier algorithm with a variable
%   step so as to have a maximum nonlinear phase rotation into each step
%   equal to X.dphimax. However, the step cannot be larger than X.dzmax. 
%
%   The CNLSE uses the same rules except that each X.length/X.nplates
%   meters a new PMF is generated.
%
%   Alternatively, the step can be chosen adaptively basing the choice on a
%   target local truncation error (NLSE only). In such a case the following 
%   parameters should be added to X:
%
%       X.ltol:       local truncation error, i.e. max distance between the
%                     field obtained by moving once or twice in a step.
%       X.dphiadapt:  true/false. True: the local truncation error method
%                     is applied only in the first step and used to correct 
%                     X.dphimax. After the first step the SSFM proceeds  
%                     using the approach based on X.dphimax. Default: false
%
%   FLAG is a string of four characters governing the type of propagation. 
%   The first character is 'g' if GVD (i.e. beta2,beta3) is on or '-' in 
%   absence of GVD. Note that with 'sepfields' in reset_all.m this function 
%   accounts for the walkoff effect even with the GVD flag set to '-'. 
%   The second character is 'p' for propagation of a polarized field in 
%   presence of birefringence and PMD or '-' in absence of such effects. 
%   The third is 's' if SPM is on or '-' in absence of SPM. 
%   Likewise, the fourth character is 'x' or '-' in presence/absence of 
%   XPM. The most complete case is FLAG='gpsx' and corresponds to 
%   propagation in presence of fiber GVD+PMD+SPM+XPM. 
%
%   E.g. Propagation with GVD+XPM, without PMD and SPM -> FLAG='g--x'
%
%   The fourth character of FLAG is active only with channels separated 
%   (see option 'sepfields' in CREATE_FIELD). In this case, the propagation 
%   neglects the effect of four-wave mixing, which can be taken in account 
%   only by combining all channels into a unique field and hence it is a 
%   special case of SPM.
%
%   OUT=FIBER(X,FLAG) returns in OUT a struct containing the birefringence
%   parameters used by fiber:
%   
%   OUT.db0 = birefringence [rad] at GSTATE.FN=0 (see RESET_ALL).
%   OUT.theta = azimuth [rad] of all the PMFs composing the fiber.
%   OUT.epsilon = ellepticity [rad] of all the PMFs composing the fiber.
%   OUT.dgd = DGD [symbols].
%   OUT.lcorr = length [m] of each PMF trunk.
%   OUT.betat = beta(omega), i.e. scalar phase shift [rad] including GVD,
%       slope,etc, where omega/2/pi is the vector of FFT frequencies. betat
%       is common to both polarizations.
%   OUT.db1 = differential phase shift [rad] induced by PMD.
%
%   OUT can be used to recover the PMD transfer matrix of the fiber (see 
%   INVERSE_PMD).
%
%   The CNLSE is described as the concatenation of X.nplates PMF trunks,
%   each with principal states of polarization randomly distributed over the 
%   Poincare sphere, if not otherwise specified. Each PMF has constant DGD 
%   and randomly distributed birefringence. The nonlinearity is inserted 
%   after a certain number of trunks, depending on X.dzmax and X.dphimax. 
%   The diagram is the following:
%
%    ------------ -----------     --------- ------------ -----------
%   |  Lin step  | Lin step  |   | NL STEP |  Lin step  | Lin step  |
%   |    PMF 1   |   PMF 2   |...|         |    PMF k   |   PMF k+1 |...
%    ------------ -----------     --------- ------------ -----------
%
%   where PMF k, k=1,2,... is a PMF fiber randomly chosen on the Poincare
%   sphere.
%
%   See also FIBER_GUI, INVERSE_PMD
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
CLIGHT = CONSTANTS.CLIGHT;      % speed of lightin vacuum [m/s]
global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

SAFETYFCT = 0.9;    % Safety step-reduction factor
DEF_PLATES = 100;   % number of waveplates when birefringence is on
[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
Nfft = GSTATE.NSYMB*GSTATE.NT;
nfnames = fieldnames(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK AND INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    error('Missing propagation type');
end
if ~any(strcmp(nfnames,'dzmax')) || (x.dzmax > x.length)
    x.dzmax = x.length; 
end

if any(strcmp(nfnames,'ltol')), % adaptive step-size
    if ~any(strcmp(nfnames,'dphimax')), x.dphimax = Inf;end
    trg.err = x.ltol;
    trg.safety = SAFETYFCT;
    if any(strcmp(nfnames,'dphiadapt')) && x.dphiadapt
        tolflag = 1; % adaptive first step, then constant phase method
    else
        tolflag = 2; % adaptive steps with target local error
    end
else
    trg = [];
    tolflag = 0; % constant phase method
end

fls = [0 0 0 0];  % flag converted into array of numbers
switch lower(flag)
    case '----'  % only attenuation
        dphimaxt = Inf;
        dzmaxt = x.length;
    case 'g---'  % only GVD
        fls(1) = 1;
        dphimaxt = Inf;
        dzmaxt = x.length;
    case '-p--'  % only PMD
        fls(2) = 1;
        dphimaxt = Inf;
        dzmaxt =x.length;        
    case '--s-'  % only SPM
        fls(3) = 1;
        if nfc == 1
            dphimaxt = Inf; % with one field we have the exact solution
            dzmaxt = x.length;
        else
            dphimaxt = x.dphimax;
            dzmaxt = x.dzmax;     
        end
    case '---x'  % only XPM 
        if nfc == 1
            error('flag ''---x'' available only for channels separated');
        end
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;              
        fls(4) = 1;
    case 'gp--'  % GVD+PMD        
        fls(1:2) = 1;
        dphimaxt = Inf;
        dzmaxt = x.length;                       
    case 'g-s-'  % GVD+SPM
        fls([1 3]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;              
    case 'g--x'  % GVD+XPM
        if nfc == 1
            error('flag ''g--x'' available only for channels separated');
        end        
        fls([1 4]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;    
    case '-ps-'  % PMD+SPM
        fls([2 3]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;              
    case '-p-x'  % PMD+XPM
        if nfc == 1
            error('flag ''-p-x'' available only for channels separated');
        end        
        fls([2 4]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;              
    case '--sx'  % SPM+XPM
        fls(3) = 1;
        if (nfc ~= 1)
            fls(4) = 1;
            dphimaxt = x.dphimax;
            dzmaxt = x.dzmax; 
        else
            dphimaxt = Inf;
            dzmaxt = x.length;
        end
    case 'g-sx'  % GVD+SPM+XPM
        fls([1 3]) = 1;
        if (nfc ~= 1), fls(4) = 1;end;        
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;     
    case '-psx'  % PMD+SPM+XPM
        fls(2:3) = 1;
        if (nfc ~= 1), fls(4) = 1;end;        
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;
    case 'gps-'  % GVD+PMD+SPM
        fls(1:3) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;     
    case 'gp-x'  % GVD+PMD+XPM
        if nfc == 1
            error('flag ''gp-x'' available only for channels separated');
        end     
        fls([1 2 4]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;             
    case 'gpsx'  % GVD+PMD+SPM+XPM
        fls(1:3) = 1;
        if (nfc ~= 1), fls(4) = 1;end;        
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;     
        
    otherwise
        error('wrong flag. E.g. flag can be ''g---'',''gp--'',''-s--'', etc');
end

isy = ~isempty(GSTATE.FIELDY); 
isv = (fls(2) == 1) || isy;
if fls(2) == 1
    if ~any(strcmp(nfnames,'manakov'))
        x.manakov = 'no';
    end 
    if ~any(strcmp(nfnames,'dgd')), error('Missing DGD in fiber');end
    isdb0 = any(strcmp(nfnames,'db0'));         % exist db0?
    istheta = any(strcmp(nfnames,'theta'));
    isepsilon = any(strcmp(nfnames,'epsilon'));
    ispmf = isdb0+istheta+isepsilon;
    if ispmf == 3   % PMF fiber
        x.nplates = length(x.theta);
        brf.db0 = x.db0;          % delta beta 0 of a PMF fiber
        brf.theta=x.theta;        % theta and epsilon are the angles of the PMF
        brf.epsilon=x.epsilon;
        dgdrms = x.dgd/x.nplates; % DGD per trunk
    elseif ispmf == 0 
        if ~any(strcmp(nfnames,'nplates'))
            x.nplates = DEF_PLATES;    % default value
        end
        brf.db0 = rand(x.nplates,1)*2*pi - pi;   % general case: No PMF.   
        brf.theta=rand(x.nplates,1)*pi - 0.5*pi;     % azimuth: uniform R.V.;   
        brf.epsilon=0.5*asin(rand(x.nplates,1)*2-1); % uniform R.V. over the Poincare sphere       
        dgdrms = sqrt((3*pi)/8)*x.dgd/sqrt(x.nplates);
        % DGD rms per trunk. x.dgd is the average value. 3pi/8 comes from
        % the Maxwellian distribution of DGD
    else                                             
        error('Missing one of db0, theta or epsilon in fiber');
    end
    brf.dgd = x.dgd;
    dgdrms = dgdrms/GSTATE.SYMBOLRATE; % convert in [ns]
    if ~isy, 
        GSTATE.FIELDY = zeros(nfr,nfc); % create y component
        warning('optilux:fiber',['You are working with two polarizations ',...
            'but in create_field you initialized just one polarization']);
    end
else
    dgdrms = 0; % turn off all polarization and birefringence effects
    brf.db0 = 0;
    brf.theta = 0;
    brf.epsilon = 0;
    ispmf = 1;
    x.manakov = 'no';
    x.nplates = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alphalin = (log(10)*1e-4)*x.alphadB;      % [m^-1]
if alphalin == 0
    Leff = x.length;      % effective length [m]
else
    Leff = (1-exp(-alphalin*x.length))/alphalin;
end
b20 = -x.lambda^2/2/pi/CLIGHT*x.disp*1e-6; % beta2 [ns^2/m] @ lambda
b30 = (x.lambda/2/pi/CLIGHT)^2*(2*x.lambda*x.disp+x.lambda^2*x.slope)*1e-6;  
                                     % beta3 [ns^3/m] @ lambda
b30 = b30*fls(1);   % insert b30 only with GVD flag on.                                     

maxl = max(GSTATE.LAMBDA);
minl = min(GSTATE.LAMBDA); 
lamc = 2*maxl*minl/(maxl+minl); %central wavelength: 1/lamc = 0.5(1/maxl+1/minl)

% Domega_ik: [1/ns]. "i" -> at ch. i, "0" -> at lambda, "c" -> at lamc
Domega_i0 = 2*pi*CLIGHT*(1./GSTATE.LAMBDA-1/x.lambda);    
Domega_ic = 2*pi*CLIGHT*(1./GSTATE.LAMBDA-1/lamc);    
Domega_c0 = 2*pi*CLIGHT*(1./lamc-1/x.lambda); 
b1 = b20*Domega_ic+0.5*b30*(Domega_i0.^2-Domega_c0.^2);  %ch's beta1 [ns/m]
if nfc == 1
    beta1 = 0;    % [ns/m] @ lamc
    Domega_i0 = 2*pi*CLIGHT*(1./lamc-1/x.lambda);        
	gam = 2*pi*x.n2./(lamc*x.aeff)*1e18; % nonlinear index [1/mW/m]
else                                     
    beta1 = b1;   % [ns/m] @ GSTATE.LAMBDA
	gam = 2*pi*x.n2./(GSTATE.LAMBDA*x.aeff)*1e18; % nonlinear index [1/mW/m]
end
beta2 = b20+b30*Domega_i0;  % beta2 [ns^2/m]@ lamc (nfc=1) or GSTATE.LAMBDA
                            % (nfc ~= 1).
beta2 = beta2*fls(1);    % insert GVD only if the GVD flag is on. 
% Note:  fls(1)=0 means that we are neglecting local GVD, but we account
% for the walk-off effect.

Dch = x.disp+x.slope*(GSTATE.LAMBDA-x.lambda);  % dispersion of the channels
% Note: channel 1 is the reference channel. All the other channels shift
%   (beta1 effect) with respect to it.

Ld = Inf*ones(1,GSTATE.NCH);
fzb2 = find(Dch);
Ld(fzb2) = 1./(GSTATE.SYMBOLRATE^2*abs(x.lambda^2/2/pi/CLIGHT*Dch(fzb2)*1e-6));
% dispersion length [m]
if b30 ~= 0
    Lds = 1./(GSTATE.SYMBOLRATE^3*abs(b30));      % slope length [m]
else
    Lds = Inf;
end

betat = zeros(Nfft,nfc);    % beta, whose Taylor coeffs are the betai
db1 = zeros(Nfft,nfc);  % delta beta, i.e. phase shift of DGD.
omega = 2*pi*GSTATE.SYMBOLRATE*GSTATE.FN';     % angular frequency [rad/ns]

for kch=1:nfc
    betat(:,kch) = omega*beta1(kch)+0.5*omega.^2*beta2(kch)+...
        omega.^3*b30/6;   % beta coefficient [1/m]
    if fls(2) == 1
        db1(:,kch) = dgdrms*omega;  % DGD phase shift  
        % since the step is non-uniform, the DGD will be normalized in
        % matrix_step
    end
end

Lnl = 1./(gam.*GSTATE.POWER);    % nonlinear length [m]

%%%%%%%%%%%%%%%%%% UPDATE DELAY AND CUMULATED DISPERSION %%%%%%%%%%%%%%%%%%
loc_delay = x.length*GSTATE.SYMBOLRATE.*b1;
GSTATE.DELAY = GSTATE.DELAY+ones(1+isy,1)*loc_delay;
GSTATE.DISP = GSTATE.DISP + ones(1+isy,1)*fls(1)*Dch*x.length*1e-3;% cum .dispersion [ps/nm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SSFM PROPAGATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%
if tolflag == 2 % adaptive step, local truncation error method
    if isv  % with polarization
        error('adaptive step available in absence of polarization effects');
    else    % scalar propagation
        [firstdz,ncycle,GSTATE.FIELDX]=scalar_a_ssfm(GSTATE.FIELDX,betat,...
            dzmaxt,dphimaxt,gam,alphalin,Nfft,nfc,x.length,trg,fls);
    end
else            % constant phase x step
    if isv  % with polarization
        [firstdz,ncycle,GSTATE.FIELDX,GSTATE.FIELDY,brf]=matrix_ssfm(...
            GSTATE.FIELDX,GSTATE.FIELDY,betat,db1,dzmaxt,dphimaxt,gam,...
            alphalin,nfc,x.length,x.nplates,x.manakov,fls,brf);
        if nargout, zbrf=brf;end
    else    % scalar propagation
        [firstdz,ncycle,GSTATE.FIELDX]=scalar_ssfm(GSTATE.FIELDX,betat,...
            dzmaxt,dphimaxt,gam,alphalin,Nfft,nfc,x.length,fls,tolflag,trg);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRINT SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GSTATE.PRINT

    if nfc == 1
        gamprint = gam*ones(1,GSTATE.NCH);
    else
        gamprint = gam;
    end

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===              fiber               ===\n');
    fprintf(fid,'========================================\n\n');
    fprintf(fid,'Fiber parameters:\n\n');
    fprintf(fid,'Length:%17.3f  [km]\n',x.length*1e-3);
    fprintf(fid,'Attenuation:%12.2f  [dB/km] (Leff = %7.3f [km])\n',...
	x.alphadB,Leff*1e-3);
    fprintf(fid,'lambda of Dc:%11.2f  [nm]\n',x.lambda);
    fprintf(fid,'Dc:%21.4f  [ps/nm/km]\n',x.disp);
    fprintf(fid,'Slope:%18.4f  [ps/nm^2/km]\n',x.slope);
    fprintf(fid,'n2:%21.2e  [m^2/W]\n',x.n2);
    fprintf(fid,'Aeff:%19.2f  [um^2]\n\n',x.aeff);
    if fls(2) == 1
        fprintf(fid,'DGD:%12.4f  [bits]\n',x.dgd);
        fprintf(fid,'# plates:%d  \n',x.nplates);
        fprintf(fid,'Manakov Equation: %s\n',x.manakov);
        if ispmf == 3
            fprintf(fid,'db0 = %8.2f, theta = %3.2f*pi, epsilon = %3.2f*pi\n',...
                brf.db0(1),brf.theta(1)/pi,brf.epsilon(1)/pi);
        else
            fprintf(fid,'Random birefringence\n');
        end
    end        
    fprintf(fid,'Propagation type: ''%s''\n\n',flag);
    if tolflag, fprintf(fid,'Local error x step: %.1e\n',x.ltol); end
    fprintf(fid,'Max NL phase rotation x step: %-6.2g  [rad]\n',dphimaxt);
    fprintf(fid,'Max step: %.2e  [m]\n',dzmaxt);
    fprintf(fid,'Initial step: %.2e  (num. steps: %d)\n\n',firstdz,ncycle);
    fprintf(fid,'Channel properties (Ld: disp. length. Lnl: NL length):\n\n');
    for kch=1:GSTATE.NCH
	fprintf(fid,'ch. #%.2d: Dc = %.4f  [ps/nm/km]   (Ld = %3.2e [km])\n',...
            kch,Dch(kch),Ld(kch)*1e-3);
	fprintf(fid,'\t gamma = %.3e [1/mW/km] (Lnl = %3.2e [km])\n',...
            gamprint(kch)*1e3,Lnl(kch)*1e-3);
	fprintf(fid,'\t sqrt(Ld/Lnl) = %.4f\n',sqrt(Ld(kch)/Lnl(kch)));
	fprintf(fid,'\t local delay = %.3f\n',loc_delay(kch));
    end
    fprintf(fid,'\nSlope length Lds: %-3.2e  [km]\n',Lds*1e-3);
    fprintf(fid,'\nGlobal  delay (ch.1 -> %d)\n',GSTATE.NCH);
    for kch=1:GSTATE.NCH
	fprintf(fid,'%.3f  ',GSTATE.DELAY(kch));
    end
    fprintf(fid,'\nGlobal cumulated dispersion [ps/nm] ');
    fprintf(fid,'(ch.1 -> %d)\n',GSTATE.NCH);
    for kch=1:GSTATE.NCH
	fprintf(fid,'%.3f  ',GSTATE.DISP(kch));
    end

    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT

%--------------------------------------------------------------------------
function [firstdz,ncycle,ux,uy,brf]=matrix_ssfm(ux,uy,betat,db1,dzmaxt,...
    dphimaxt,gam,alphalin,nfc,Lf,nplates,manakov,fls,brf)

%MATRIX_SSFM SSFM algorithm for the CNLSE with PMD
%   [FIRSTDZ,NCYCLE,UX,UY,BRF]=MATRIX_SSFM(UX,UY,BETAT,DB1,DZMAXT,...
%   DPHIMAXT,GAM,ALPHALIN,NFC,LF,NPLATES,MANAKOV,FLS,BRF) solves the Coupled
%   Nonlinear Schroedinger equation through the split step Fourier method
%   (SSFM). 
%
%   UX and UY are the x and y components of the electric field,
%   respectively.
%   BETAT is the the (scalar) beta(omega) coefficient, DB1 is the 
%   differential phase shift induced by PMD, with DB1=DB1(omega). 
%   DZMAXT is the largest step allowed by nonlinear effects. DPHIMAXT is 
%   the largest nonlinear phase rotation x step. GAM is the nonlinear 
%   coefficient [1/mW/m].
%   ALPHALIN is the attenuation [m^-1]. NFC is the number of fields. LF is
%   the fiber length [m]. NPLATES is the number of PMD emulation
%   waveplates. MANAKOV='yes' or 'no' for solving the Manakov equation
%   (faster) instead of the CNLSE. FLS is the propagation flag.
%   
%   BRF.db0, BRF.theta, BRF.epsilon are the birefringence at frequency=0, 
%   the azimuth and the ellepticity, respectively. They are vectors of 
%   length equal to NPLATES.
%
%   In every step the fiber is supposed to be a PMF with principal states 
%   of polarization (PSPs) oriented as indicated by BRF.theta and 
%   BRF.epsilon. The DGD is assumed to be independent from omega. 
%
%   The function returns UX and UY after propagation, the first step 
%   FIRSTD, the number of cycles NCYCLE used by the SSFM, and adds new
%   fields to BRF.

% Pauli matrices
sig0 = eye(2);
sig2 = [0 1;1 0];
sig3i = [0 1;-1 0]; % =i*sig3 = i*[0 -i;i 0]
%%%%%%%%%%%%%%%

Nfft = size(ux,1);
if strcmp(manakov,'yes') 
    gam=gam*8/9;    % 8/9 factor of Manakov equation
    ismanakov = true;
else
    ismanakov = false;
end % Manakov equation

ncycle = 1;      % number of cycles
lcorr = Lf/nplates; % waveplate length [m] (PMD step)
dz_miss = 0;    % non-propagated distance within the waveplate from previous
                % run.
brf.lcorr = lcorr;

dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,ux,uy,true);

halfalpha = 0.5*alphalin; 
ntot = 0; % In which waveplate am I now?
firstdz = dz;
zprop = dz;                   % running distance [m]
while zprop < Lf        % all steps except the last
    
    [ux,uy] = matrix_nl_step(Nfft,ismanakov,...
        alphalin,gam,dz,ux,uy,nfc,fls(3),fls(4));

    %%% GVD+PMD %%%
    [dzb,dz_miss,nmem,ntrunk] = checkstep(zprop,dz,lcorr,dz_miss,ntot);    
    
    [ux,uy] = matrix_step(betat,db1,dzb,ntrunk,ux,uy,...
        sig0,sig2,sig3i,brf,ntot,nmem);
   
    ntot = ntot + ntrunk - nmem; % update waveplate number
      
    ux = ux*exp(-halfalpha*dz);    %  attenuation
    uy = uy*exp(-halfalpha*dz);    %  attenuation
    
    dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,ux,uy,true);
    zprop = zprop+dz;
    ncycle = ncycle+1;
end 
last_step = Lf-zprop+dz; % last step
%    fprintf('last dz=%8.2f dzbnext=%8.2f nstep_nl=%5d\n',last_step,dzb,nstep_nl) 

[ux,uy] = matrix_nl_step(Nfft,ismanakov,alphalin,gam,...
    last_step,ux,uy,nfc,fls(3),fls(4)); % last NL
                                                          
%%% GVD+PMD %%%
[dzb,dz_miss,nmem,ntrunk] = checkstep(Lf,last_step,lcorr,dz_miss,ntot);

[ux,uy] = matrix_step(betat,db1,dzb,ntrunk,ux,uy,sig0,sig2,sig3i,...
    brf,ntot,nmem); % last LIN

ux = ux*exp(-halfalpha*last_step);    % attenuation
uy = uy*exp(-halfalpha*last_step);    % attenuation

brf.betat = betat;
brf.db1 = db1;

%--------------------------------------------------------------------------
function [firstdz,ncycle,u]=scalar_ssfm(u,betat,dzmaxt,dphimaxt,gam,...
    alphalin,Nfft,nfc,Lf,fls,tolflag,trg)

%SCALAR_SSFM SSFM algorithm for the scalar NLSE 
%   [FIRSTDZ,NCYCLE,U]=SCALAR_SSFM(U,BETAT,DZMAXT,DPHIMAXT,GAM,ALPHALIN,...
%   NFFT,NFC,LF,FLS,TOLFLAG,TRG)  
%   solves the Nonlinear Schrodinger equation through the split step Fourier 
%   method (SSFM). This function works for a scalar field U.
%
%   BETAT is the the (scalar) beta(omega) coefficent. DZMAXT is the
%   largest step allowed by nonlinear effects. DPHIMAXT is the largest
%   nonlinear phase rotation x step. GAM is the nonlinear coefficient.
%   ALPHALIN is the attenuation [m^-1]. NFFT is the number of FFT points.
%   NFC is the number of fields. LF is the fiber length [m]. 
%   FLS is the propagation flag. TOLFLAG=1 corrects DPHIMAXT adaptively.
%   TRG is a struct containing the flags of the adaptive step algorithm.
%
%   This function operates over the global field GSTATE.FIELDX in absence 
%   of PMD effects. In presence of PMD see MATRIX_SSFM.
%
%   The function uses the symmetric SSFM if the step is chosen adaptively 
%   on the basis of the local truncation error. Otherwise a basic SSFM is
%   used.
%
%   The function returns the electric field U after propagation, the
%   first step FIRSTDZ and the number of cycles NCYCLE used in the SSFM.

gamrep = repmat(gam,Nfft,1);
dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,u,[],false);
halfalpha = 0.5*alphalin; 
ncycle = 1;      % number of cycles
if tolflag == 1 % adaptive search of dphimax
    if dz >= dzmaxt
        Umax = max(real(u).^2+imag(u).^2);
        maxpow = max(gam.*Umax);
        if alphalin == 0
            dphimaxt = maxpow*dz;
        else
            dphimaxt = maxpow*(1-exp(-alphalin*dz))/alphalin;
        end
    end
    dzini = dz;
    zdone=0;
    while zdone == 0
        [u ncycle nrej zdone dz]=adaptssfm(u,zdone,dz,alphalin,...
            gamrep,nfc,fls,betat,halfalpha,trg,0,0);
    end
    if dz > dzmaxt, dz = dzmaxt;end
    % on the basis of the first step correct the max phase rotation x step
%     fprintf('dphiINI=%.5f\n',dphimaxt)
    dphimaxt = dphimaxt*(1-exp(-alphalin*zdone))/(1-exp(-alphalin*dzini));
%      fprintf('dphiADAPT= %.5f\n',dphimaxt); 
    firstdz = zdone;
    zprop = zdone+dz;
    ncycle = ncycle+1;
else
    firstdz = dz;
    zprop = dz;             % running distance [m]
end
while zprop < Lf        % all steps except the last
    u = nl_step(alphalin,gamrep,dz,u,nfc,fls(3),fls(4));
                                                     % 1/3) NON-linear step
    
    u = lin_step(betat*dz,u); % 2/3) linear step    
    
    u = u*exp(-halfalpha*dz);    % 3/3) attenuation
    
    dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,u,[],false);
    zprop = zprop+dz;
%     semilogy(ncycle,dz,'b+')
    ncycle = ncycle+1;
%     fprintf('zdone = %.3f\t  dz=%.3f\t \t  ncycle=%d\n',zprop,dz,ncycle)

end 
last_step = Lf-zprop+dz;
u = nl_step(alphalin,gamrep,last_step,u,nfc,fls(3),fls(4)); 
                                                          % NON-linear step
u = lin_step(betat*last_step,u); % linear step    

u = u*exp(-halfalpha*last_step);    % attenuation

%--------------------------------------------------------------------------
function [firstdz,ncycle,u]=scalar_a_ssfm(u,betat,dzmaxt,dphimaxt,gam,...
    alphalin,Nfft,nfc,Lf,trg,fls)

%SCALAR_A_SSFM SSFM algorithm for the scalar NLSE with adaptive step-size 
%   [FIRSTDZ,NCYCLE,U]=SCALAR_A_SSFM(U,BETAT,DZMAXT,DPHIMAXT,GAM,...
%       ALPHALIN,NFC,LF,TRG,FLS)  
%   solves the Nonlinear Schrodinger equation through the split step Fourier 
%   method (SSFM) using adaptive step-size. This function works for a 
%   scalar field U.
%
%   BETAT is the the (scalar) beta(omega) coefficent. DZMAXT is the
%   largest step allowed by nonlinear effects. DPHIMAXT is the largest
%   nonlinear phase rotation x step. GAM is the nonlinear coefficient.
%   ALPHALIN is the attenuation [m^-1]. Nfft is the number of FFT points.
%   NFC is the number of fields. LF is the fiber length [m]. TRG is a 
%   struct containing in TRG.ERR the local truncation error and in 
%   TRG.SAFETY a safety reduction facto for the next step proposal.
%   FLS is the propagation flag, see fiber.m
%
%   This function operates over the global field GSTATE.FIELDX in absence 
%   of PMD effects. 
%
%   The function returns the electric field U after propagation, the
%   first step FIRSTDZ and the number of cycles NCYCLE used in the SSFM.

gamrep = repmat(gam,Nfft,1);
ncycle = 1;      % number of cycles
dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,u,[],false); % first step
halfalpha = 0.5*alphalin; 
firstdz = dz;
zdone = 0;
nrej = 0;
while zdone < Lf        
    if zdone + dz > Lf
        dz = Lf - zdone; % last step
    end
    [u ncycle nrej zdone dz]=adaptssfm(u,zdone,dz,alphalin,...
    gamrep,nfc,fls,betat,halfalpha,trg,nrej,ncycle);
    if dz > dzmaxt, dz = dzmaxt;end

end 

%--------------------------------------------------------------------------
function dz_nl=nextstep(dzmax,phimax,gam,alphalin,ux,uy,isv)

%NEXTSTEP step for the SSFM algorithm
%   DZ=NEXTSTEP(DZMAX,PHIMAX,GAM,ALPHALIN,U,NFC,ISV) evaluates the step of  
%   the SSFM algorithm for a fiber having nonlinear gamma coefficient GAM 
%   [1/mW/km] and attenuation ALPHALIN [m^-1]. U is the matrix of electrical
%   fields. NFC is the number of channels (equal to 1 if all channels are 
%   combined into a unique field). ISV=1 if the y component is on.
%   The step corresponds to a maximum nonlinear phase rotation equal to 
%   PHIMAX, under the constraint that it cannot be greater than DZMAX.

if isv % Umax is the largest normalized power
    Umax = max(real(ux).^2+imag(ux).^2+real(uy).^2+imag(uy).^2); 
else
    Umax = max(real(ux).^2+imag(ux).^2); 
end    
Pmax = max(gam.*Umax);   % largest gamma*power
leff = phimax/Pmax;             % effective length of the step
dl = alphalin*leff;             % ratio effective length/attenuation length

if dl >= 1
    dz_nl = dzmax;
else
    if alphalin == 0
        step = leff;
    else
        step = -1/alphalin*log(1-dl);
    end
    if step > dzmax
        dz_nl = dzmax;
    else
        dz_nl = step;
    end
end

%--------------------------------------------------------------------------
function [dzb,dz_miss,nmem,ntrunk] = checkstep(zprop,dz,lcorr,dz_miss,nz_old)

%CHECKSTEP evaluates the birefringence step
%   [DZB,DZ_MISS,NMEM,NTRUNK]=CHECKSTEP(ZPROP,DZ,LCORR,DZ_MISS,NZ_OLD) 
%   given the nonlinear step DZ evaluates the birefringence steps DZB 
%   inside DZ. 
%   LCORR is the waveplate (trunk) length, i.e. the length over which the 
%   birefringence is supposed constant and the fiber a PMF. 
%   DZ_MISS is the missing distance from the actual coordinate to the end 
%   of the actual birefringence waveplate. On output DZ_MISS is updated 
%   after propagating of DZ. NZ_OLD is the trunk number before applying the 
%   step.
%
%   On output DZB is a vector containing all the steps over which the
%   linear birefringence operates. NTRUNK is the number of birefringence 
%   steps within DZ. Hence length(DZB) = NTRUNK. NMEM=0 means that DZB(1)
%   starts with a new trunk, while NMEM=1 means that DZB(1) starts in the
%   same trunk as the previous step
%
%   All lengths are in [m]. 

nz = zprop/lcorr;  
nzc = ceil(nz);     % trunk number after propagation
if dz_miss == 0     % are you at the end of the trunk?
   nmem = 0;        % do not keep memory of the trunk seed
   ntrunk = nzc - nz_old;   % number of trunks within dz
   dzlast = dz-lcorr*(ntrunk-1);
   dzb = [lcorr*ones(1,ntrunk-1),dzlast];
   dz_miss = lcorr - dzlast; % update missing coordinate within the trunk
else
   nmem = 1;        % do not keep memory
   ntrunk = nzc - nz_old + 1;
   if ntrunk == 1
       dzb = dz;
       dz_miss = dz_miss-dz;
   else % dz moves across >1 trunks
       dzlast = dz-dz_miss-lcorr*(ntrunk-2);
       dzb = [dz_miss,lcorr*ones(1,ntrunk-2),dzlast];
       dz_miss = lcorr - dzlast;
   end    
end 


%-----------------------------------------------------------
function u = lin_step(betaxdz,u)

%LIN_STEP linear fiber in scalar propagation
%   U=LIN_STEP(BETAXDZ,U) propagates the electric field U into a purely
%   linear fiber having beta*dz factor BETAXDZ, i.e. a fiber
%   with transfer function exp(-i*beta(omega)*dz), where  betat is
%   approximated by its Taylor series in omega.


Hf = fastexp(-betaxdz); %  Fast exponential: fastexp(x) = exp(i*x)

u = ifft( fft(u) .* Hf);

%-----------------------------------------------------------
function u = nl_step(alphalin,gam,dz,u,nfc,spm,xpm)

%NL_STEP nonlinear fiber in scalar propagation
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
pow = (real(u).^2+imag(u).^2);
if xpm
    if spm
        pow = 2*sum(pow,2)*ones(1,nfc)-pow; % sum(.,2): row sum
    else
        pow = 2*(sum(pow,2)*ones(1,nfc)-pow);        
    end
else
    if ~spm
        return
    end
end
u = u.*fastexp(-gam.*pow*leff);

%-----------------------------------------------------------
function [ux,uy] = matrix_nl_step(Nfft,ismanakov,alphalin,gam,dz,ux,uy,...
    nfc,spm,xpm)

%MATRIX_NL_STEP nonlinear fiber of the CNLSE
%   [UX,UY] = MATRIX_NL_STEP(NFFT,ISMANAKOV,ALPHALIN,GAM,DZ,UX,UY,NFC,SPM,XPM) 
%   propagates the electric field of x and y components UX and UY,
%   respectively, into a purely nonlinear fiber having nonlinear 
%   coefficient GAM [1/mW/km], attenuation ALPHA_LIN [m^-1]. NFC 
%   is the number of channels. DZ is the fiber length [m]. SPM and XPM are 
%   1 or 0 if self-phase modulation and cross-phase modulation have to be 
%   taken in account, respectively. NFFT is the number of FFT points x
%   channel. ISMANAKOV=1 if the nonlinear step is the one of the Manakov 
%   equation. For spm only this function implements the following:
%
%   U = expm(-i*gam*leff*(|U|^2*[1 0;0 1] - 1/3*U'sig3*U*sig3))
%
%   where U = [ux(k);uy(k)], k is the kth frequency, k=1,2,...,Nfft, and 
%   sig3 = [0 -i;i 0]. leff is the fiber effective length [m].


if alphalin == 0
    leff = dz;
else
    leff = (1-exp(-alphalin*dz))/alphalin;  % effective length of dz
end
for kch=1:nfc
    if spm          %  Fast exponential: fastexp(x) = exp(i*x)
        power = real(ux(:,kch)).^2+imag(ux(:,kch)).^2+...
            real(uy(:,kch)).^2+imag(uy(:,kch)).^2;
        gamleff = gam(kch)*leff;
        nlscalar = fastexp(-gamleff*power);
        ux(:,kch) = ux(:,kch) .* nlscalar;
        uy(:,kch) = uy(:,kch) .* nlscalar;
        % 2: -> real part. 1/3: weight of non linear term in CNLSE
        if ~ismanakov % solve the CNLSE
            s3 = 2*(real(ux(:,kch)).*imag(uy(:,kch))-...
                imag(ux(:,kch)).*real(uy(:,kch))); % stokes comp. #3
            expigamleffs3d3 = fastexp(gamleff*s3/3);
            cosphi = real(expigamleffs3d3); % (fast) cos(gamleff*s3/3) 
            sinphi = imag(expigamleffs3d3); % (fast) sin(gamleff*s3/3)
            uxx = cosphi.*ux(:,kch) + sinphi.*uy(:,kch);
            uyy = -sinphi.*ux(:,kch) + cosphi.*uy(:,kch);
            ux(:,kch) = uxx;
            uy(:,kch) = uyy;
        end
    end
    if xpm
        error('The CNLSE with separate fields is not yet implemented');
        n = [1:kch-1,kch+1:nfc];    % XPM channels
        s3 = zeros(Nfft,1);
        for nch=1:nfc-1
            power = real(ux(:,n(nch))).^2+imag(ux(:,n(nch))).^2+...
                real(uy(:,n(nch))).^2+imag(uy(:,n(nch))).^2;
            s3 = 4/3*(real(ux(:,n(nch))).*imag(uy(:,n(nch)))-...
                imag(ux(:,n(nch))).*real(uy(:,n(nch)))) + s3; % prop. to stokes comp. #3
        % 4=2*2: -> real part * xpm factor. 1/3: weight of non linear term in CNLSE
            nlscalar = fastexp(-2*gam(n(nch))*power*leff); % *2: XPM
            ux(:,kch) = ux(:,kch) .* nlscalar;
            uy(:,kch) = uy(:,kch) .* nlscalar;            
        end
        if ismanakov
            uxx = cos(leff*s3).*ux(:,n(nch)) - sin(leff*s3).*uy(:,n(nch));
            uyy = sin(leff*s3).*ux(:,n(nch)) + cos(leff*s3).*uy(:,n(nch));
            ux(:,kch) = uxx;
            uy(:,kch) = uyy;
        end
    end
end

%--------------------------------------------------------------------------
function [ux,uy]=matrix_step(betat,db1,dzb,ntrunk,ux,uy,...
    sig0,sig2,sig3i,brf,ntot,nmem)

%MATRIX_STEP Matrix linear step
%   [UX,UY]=MATRIX_STEP(BETAT,DB1,DZB,NTRUNK,UX,UY,SIG0,SIG2,SIG3I,BRF,NTOT,NMEM)  
%   applies the linear vectorial step including birefringence and PMD. 
%   BETAT is the scalar part of the beta(omega) coefficient, i.e. 
%   BETAT= beta0 + beta1*omega + beta2/2*omega^2+beta3/6*omega^3
%
%   DZB is the birefringence step [m], while NTRUNK are the number of 
%   DZBs inside each fiber waveplate. 
%
%   SIG0 = [1 0;0 1], SIG2 = [0 1;1 0], SIG3I = [0 1;-1 0] are the Pauli 
%   matrices except SIG3I that is i times the third Pauli matrix 
%   (engineering notation).
%
%   UX and UY are the X and Y components of the electric field,
%   respectively.
%
%   BRF.db0, BRF.theta, BRF.epsilon are the birefringence at frequency=0, 
%   the azimuth and the ellepticity, respectively. DB1 is the differential 
%   phase shift induced by PMD, with DB1=DB1(omega). Hence, on the PSPs
%   basis PMD acts as exp(i*DB1) and exp(-i*DB1).
%   NTOT is the waveplate number before applying the matrix step. NMEM=1 is
%   the propagation starts in the previous waveplate, or NMEM=0 if the
%   propagation starts in a new waveplate.

ux = fft(ux);
uy = fft(uy);

for k=1:ntrunk
    n = ntot+k-nmem; % trunk number
    
    matRth = cos(brf.theta(n))*sig0 - sin(brf.theta(n))*sig3i;    % orthogonal matrix
    matRepsilon = complex(cos(brf.epsilon(n))*sig0,sin(brf.epsilon(n))*sig2); % orthogonal
    matR = matRth*matRepsilon;  % matrix of change of basis over the PSPs.

    % Note: Calling A=[GSTATE.FIELDX(k);GSTATE.FIELDY(k)] the electric 
    %   field for the kth frequency, we have that matR*D*matR'*A
    %   is the linear PMD step, where D is the diagonal matrix
    %   where the DGD operates.

    % 1) move onto the PSPs basis
    uux = conj(matR(1,1))*ux + conj(matR(2,1))*uy;
    uuy = conj(matR(1,2))*ux + conj(matR(2,2))*uy;  
    
    % 2) apply birefringence, DGD and GVD: all in a diagonal matrix    
    combeta=betat*dzb(k);           % common beta factor
    deltabeta=0.5*(db1+brf.db0(n))*dzb(k)/brf.lcorr;  % differential beta factor
    % Note: dzb(k)/brf.lcorr: fraction of DGD within current step dzb(k).
    uux = fastexp(-(combeta+deltabeta)).*uux; 
    uuy = fastexp(-(combeta-deltabeta)).*uuy;   
    
    % 3) come back in the original basis
    ux = matR(1,1)*uux + matR(1,2)*uuy;
    uy = matR(2,1)*uux + matR(2,2)*uuy;
end
ux = ifft(ux);
uy = ifft(uy);

%--------------------------------------------------------------------------
function [u ncycle nrej zdone dz]=adaptssfm(u,zdone,dz,alphalin,...
    gamrep,nfc,fls,betat,halfalpha,trg,nrej,ncycle)

% ADAPTSSFM single SSFM adaptive step in the scalar case
%   [U NCYCLE NREJ ZDONE DZ]=ADAPTSSFM(U,ZDONE,DZ,ALPHALIN,...
%       GAMREP,NFC,FLS,BETAT,HALFALPHA,TRG,NREJ,NCYCLE)
%   evaluates a single step of the SSFM algorithm using adaptive step size
%   measure. U is the electric field, NCYCLE the cumulated number of
%   cycles, NREJ the number of steps rejected, ZDONE the cumulated distance,
%   DZ the proposed step on input and the proposed next step on output.
%   ALPHALIN is the fiber attenuation in [m^-1], GAMREP the fiber gamma
%   index [1/mW/km], NFC the number of electric fields, FLS a propagation
%   flag, BETAT is the the (scalar) beta(omega) coefficent, HALFALPHA is 
%   ALPHALIN/2, TRG is a struct containing in TRG.ERR the local truncation 
%   error and in TRG.SAFETY a safety reduction facto for the next step
%   proposal.
%
%   The step is adaptively chosen using the standard technique based on
%   Richardson extrapolation. See [1]. The main difference with [1] is that
%   the local error is an absolute error instead of relative and that the
%   update rule for the step is different and based on TRG.SAFETYFCT.
%
%   [1] O. V. Sinkin, R. Holzlohner, J. Zweck and C. R. Menyuk, "
%   Optimization of the Split-Step Fourier Method in Modeling Optical-Fiber
%   Communications Systems," J. Lightw. Technol, vol. 21, n.1, pp.61-68,
%   Jan. 2003


dz2 = 0.5*dz;   % dz/2
dz4 = 0.25*dz;  % dz/4

ustack = u;     % keep memory in case of failure
uh = u;
%%% one big step (symmetric SSFM)
u = nl_step(alphalin,gamrep,dz2,u,nfc,fls(3),fls(4));
u = u*exp(-halfalpha*dz2);    % 3/3) attenuation

% 1/3) NON-linear step

u = lin_step(betat*dz,u); % 2/3) linear step
u = nl_step(alphalin,gamrep,dz2,u,nfc,fls(3),fls(4));
u = u*exp(-halfalpha*dz2);    % 3/3) attenuation


%%% two small steps
uh = nl_step(alphalin,gamrep,dz4,uh,nfc,fls(3),fls(4));
uh = uh*exp(-halfalpha*dz4);    % 3/3) attenuation
uh = lin_step(betat*dz2,uh);
uh = nl_step(alphalin,gamrep,dz2,uh,nfc,fls(3),fls(4));
uh = uh*exp(-halfalpha*dz2);    % two NL-neighboring steps

% uh = nl_step(alphalin,gamrep,dz4,uh,nfc,fls(3),fls(4));
% uh = uh*exp(-halfalpha*dz4);    % 3/3) attenuation
uh = lin_step(betat*dz2,uh);
uh = nl_step(alphalin,gamrep,dz4,uh,nfc,fls(3),fls(4));
uh = uh*exp(-halfalpha*dz4);    % 3/3) attenuation

% est_err = norm(u-uh)/dz;
% est_err = max(abs(u-uh))/dz;
est_err = max(max(sqrt(real(u-uh).^2+imag(u-uh).^2)))/dz;

if est_err > trg.err % reject the step
    dz = trg.safety*sqrt(trg.err/est_err)*dz; % reduce the step
    u = ustack;
    nrej = nrej+1;
else % accept the step
    u = 4/3*uh-1/3*u; % Richardson extrapolation
    zdone = zdone+dz;
%          semilogy(ncycle,dz,'ro')
dz = trg.safety*sqrt(trg.err/est_err)*dz; % increase the step
   ncycle = ncycle+1;
end
% fprintf('zdone = %.3f\t  dz=%.3f\t est_err=%.3e\t nrej=%d\t ncycle=%d\n',zdone,dz,est_err,nrej,ncycle)
