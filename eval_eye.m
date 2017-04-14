function varargout=eval_eye(ich,x,pat)

%EVAL_EYE Evaluate the eye opening for a non-coherent transmission.
%   VARARGOUT=EVAL_EYE(ICH,X,PAT) evaluates the eye opening of channel 
%   ICH of a non-coherent transmission. The eye opening [mA] is defined as:
%
%           max(worst1-worst0) 
%
%   being worst1 and worst0 the worst mark/space samples. The best sampling 
%   time and the max(eye) are evaluated by parabolic interpolation.
%   PAT is the symbol pattern, after decoding (see PAT_DECODER). PAT can be
%   a vector or a matrix, depending on the modulation format (see next).
%
%   X is a structure whose fields are:
%
%   X.rec = receiver type. Valid arguments are: 'ook','psbt','nf-dpsk' to 
%           use RECEIVER_OOK, 'dpsk' to use RECEIVER_DPSK, 'dqpsk' to use 
%           with RECEIVER_DQPSK.
%   X.oftype = optical filter type (see MYFILTER)
%   X.obw = optical filter 3 dB bandwidth normalized to the symbol rate
%   X.oord = optical filter order (if X.oftype is 'supergaussian')
%   X.eftype = electrical filter type (see MYFILTER)
%   X.ebw = electrical filter 3-dB bandwidth normalized to the symbol rate
%   X.eord = electrical filter order (if X.eftype is 'supergaussian')
%
%   X can also have the optional parameters:
%
%   X.ts = Fixed sampling time (-0.5 <= X.ts <= 0.5). 
%   X.plot = 'ploteye': plots the eye in the active figure; 'plotcur' plots
%       the received current.
%   X.color = color string for the plot (see PLOT). E.g. 'b-'.
%   X.dpost = post compensating fiber cumulated dispersion [ps/nm]
%   X.slopez = post compensating fiber cumulated slope  [ps/nm^2]
%   X.lambda = wavelength [nm] at which the post compensating fiber
%              has a cumulated dispersion equal to X.dpost.
%   X.comp = component on which evaluate eye and calculate BER (dqpsk 
%            modulation only). Can be 'phase' or 'quadrature' or 'both'. In
%            the last case the function gets two measurements over the
%            in-phase and quadrature components, sampled with the same
%            clock time. X.comp='both' requires PAT to be a two-column
%            matrix with the phase/quadrature binary patterns on column 1,2,
%            respectively.
%   X.print = structure for print. E.g. X.print = {'nomefile','eye'} or
%            X.print = {'nomefile','current'}, prints to file nomefile the 
%            eye or the current, respectively. nomefile will be place into
%            GSTATE.DIR within a directory ending with '.MOD'.
%   X.delay = 'theory' means that the delay uses the theoretical delay 
%           saved within GSTATE.DELAY (see CREATE_FIELD). By default the 
%           delay is measured by a cross-correlation measurement between 
%           the received current and an artificial pulse amplitude 
%           modulation (PAM) signal with ideal non-return to zero bits with
%           symbols equal to PAT. The correlation method is useful in 
%           presence of polarization mode dispersion (PMD).
%
%   The receiver is composed of an ideal, purely linear, post compensating 
%   fiber + optical filter + optical to electrical converter + electrical 
%   lowpass filter. For example see RECEIVER_OOK or RECEIVER_DPSK.
%
%   VARARGOUT can contain a variable number of argument, from 0 to 5. In 
%   the complete case it is
%
%           VARARGOUT = [EO, TS, IS, SN, DELAY]
%
%   where EO is the normalized eye opening, while TS is the best sampling 
%   time normalized to the bit time. IS is a column vector containing the 
%   sampled bits. SN is a vector containing the FFT coefficients of the 
%   signal after the optical filter (SN is used by BER_KL). DELAY is the 
%   overall system delay in bits.
%   With X.comp='both' all variables of VARARGOUT are doubled in size for
%   accounting for the in-phase and quadrature component.
%
%   It is -0.5 <= TS <= 0.5, with usually TS ~= 0.
%
%   Note 1: This function works over a copy of the electric field. All  
%           fields of the global variable GSTATE are left unchanged.
%   Note 2: All fields of X must be lowercase.
%   Note 3: Please, beware that EVAL_EYE needs the decoded pattern, not the
%   transmitted one. See the examples.
%
%   See also PATTERN, MYFILTER, BEST_EYE, RECEIVER_OOK, RECEIVER_DPSK, 
%            RECEIVER_DQPSK, BER_KL, PAT_DECODER, CORRDELAY
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


global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

ylab=''; % for plotting purposes only
switch x.rec
    case 'ook'        
        [Iric,retx] = receiver_ook(ich,x); % get the current
        [eyeb,best_ts,Iricmat,retx,delay,xopt]=geteyeinfo(ich,Iric,x,retx,pat(:));
    case 'dpsk'
        [Iric,retx] = receiver_dpsk(ich,x); % get the current
        [eyeb,best_ts,Iricmat,retx,delay,xopt]=geteyeinfo(ich,Iric,x,retx,pat(:),ylab);        
    case 'dqpsk'
        [Iric,retx] = receiver_dqpsk(ich,x); % get the current
        isfcomp = isfield(x,'comp');
        if isfcomp && strcmpi(x.comp,'both') % evalaute the eye two times   
            x.comp = 'phase';
            if isfield(x,'plot')
                subplot(1,2,1);
                ylab = ' (in-phase)';
            end
            [eyeb(1),best_ts(1),Iricmat1,retmp,delay(1),xopta(1)]=...
                geteyeinfo(ich,Iric(:,1),x,retx,pat(:,1),ylab);
            x.comp = 'quadrature';
            if isfield(x,'plot')
                subplot(1,2,2);
                ylab = ' (quadrature)';            
            end            
            [eyeb(2),best_ts(2),Iricmat2,retmp,delay(2),xopta(2)]=...
                geteyeinfo(ich,Iric(:,2),x,retx,pat(:,2),ylab);
            Iricmat=[Iricmat1;Iricmat2];
            x.comp = 'both';
            xopt=(xopta(1)+xopta(2))*0.5;   % synchronous clock equal to the mean
        elseif isfcomp && strcmpi(x.comp,'phase')
            [eyeb,best_ts,Iricmat,retx,delay,xopt]=...
                geteyeinfo(ich,Iric(:,1),x,retx,pat(:),ylab);           
        elseif isfcomp && strcmpi(x.comp,'quadrature')
            [eyeb,best_ts,Iricmat,retx,delay,xopt]=...
                geteyeinfo(ich,Iric(:,2),x,retx,pat(:),ylab);
        else
            error('X.comp must be one of ''phase'',''quadrature'',''comp''');
        end
    otherwise
        error('Flag X.mod must be ''ook'', ''dpsk'',''dqpsk''');
end

    
%%%%%%%%%%%%% output variables
if eyeb < 0, eyeb=NaN;end
if nargout >= 1, varargout(1)={eyeb};end;
if nargout >= 2, varargout(2)={best_ts};end;
if nargout >= 3
    Iric_resamp = interp1(1:GSTATE.NT,Iricmat',xopt);
    if (size(Iricmat,1) == GSTATE.NSYMB*2)
        Iric_resamp = reshape(Iric_resamp',size(Iric_resamp',1)/2,size(Iric,2))';
    end
    varargout(3)={Iric_resamp'}; % interpolate at the best sampling times
end
if nargout >= 4, varargout(4)={retx};end;
if nargout == 5, varargout(5)={delay};end;


%--------------------------------------------------------------------------
function [eyeb,best_ts,Iricmat,retx,delay,xopt]=geteyeinfo(ich,Iric,x,retx,pat,ylab)

% Eye evaluations. On input:
%   ICH: channel number
%   IRIC: electric current returned by the receiver
%   X: struct variable
%   RETX: struct variable returned by the receiver
%   PAT: symbols pattern
%   YLAB: Y-axis label for plots
%
%   On output:
%
%   EYEB: Eye opening
%   BEST_TS: best normalized sampling time (center of bit is 0)
%   RETX: struct variable used by BER_KL
%   DELAY: overall normalized (to the bit time) delay
%   XOPT: Best sampling time in discrete points.
%

global GSTATE

ipat = (pat == 1);
isy = ~isempty(GSTATE.FIELDY);

%%%% MEASURE THE DELAY

if isfield(x,'delay') && strcmp(x.delay,'theory') % theoretical delay
    if isfield(x,'b2b')
        if strcmp(x.b2b,'b2b'), avgdelay = 0;end   % back-to-back: no delay    
    else
        if (isy) 
            avgdelay = 0.5*(GSTATE.DELAY(1,ich)+GSTATE.DELAY(2,ich));
        else
            avgdelay = GSTATE.DELAY(1,ich);
        end   
    end
    delay = avgdelay+evaldelay(x.oftype,x.obw*0.5)+...
        evaldelay(x.eftype,x.ebw)+retx.post_delay;
    
else % delay estimation by cross-correlation measurement
    delay = corrdelay(Iric,ipat,GSTATE.NT,GSTATE.NSYMB);
end

halfbit = GSTATE.NT/2;

nshift = round(halfbit-delay*GSTATE.NT); % the first bit is centered at index 1
Iricmat = reshape(fastshift(Iric,nshift),GSTATE.NT,GSTATE.NSYMB*size(Iric,2))';  % Note the transpose!
min1 = min(Iricmat(ipat==1,:),[],1);  % worst mark points
max0 = max(Iricmat(ipat~=1,:),[],1);  % worst space points

if (~isempty(min1) && ~isempty(max0))
    eyeop = min1-max0;    % eye opening
    eyeop = eyeop/GSTATE.POWER(ich);
    if isfield(x,'ts') % fixed sampling time
        xopt = round((x.ts+0.5)*GSTATE.NT);
        eyeb = eyeop(xopt);
        best_ts = x.ts;
    else
        [eyeb,best_tsn] = max(eyeop);

        if GSTATE.NT == 2
            best_ts = best_tsn/GSTATE.NT-0.5;
        else    % Interpolation for the best eye opening
            if best_tsn == 1
                to = 1:3;
            elseif best_tsn == GSTATE.NT
                to = best_tsn-2:best_tsn;
            else
                to = best_tsn-1:best_tsn+1;
            end
            % Parabolic interpolation around the max using three points.
            xopt = 0.5*(to(3)^2*(eyeop(to(1))-eyeop(to(2)))+to(1)^2*(eyeop(to(2))-eyeop(to(3)))+...
                to(2)^2*(eyeop(to(3))-eyeop(to(1))))/(to(3)*(eyeop(to(1))-eyeop(to(2)))+...
                to(1)*(eyeop(to(2))-eyeop(to(3)))+to(2)*(eyeop(to(3))-eyeop(to(1))));
            eyeb = (xopt-to(2))*(xopt-to(3))/((to(1)-to(2))*(to(1)-to(3)))*eyeop(to(1))+...
                (xopt-to(1))*(xopt-to(3))/((to(2)-to(1))*(to(2)-to(3)))*eyeop(to(2))+...
                (xopt-to(1))*(xopt-to(2))/((to(3)-to(1))*(to(3)-to(2)))*eyeop(to(3));
            best_ts = xopt/GSTATE.NT-0.5;
        end
    end
else     % if (~isempty(min1) && ~isempty(max0))
    eyeb = Inf;     % the eye does not exist
    best_ts = 0;    % any time is good
    xopt = 0.5*GSTATE.NT;
end

%%%%% GRAPHICAL OUTPUT

if isfield(x,'plot')
    if ~isfield(x,'color') 
        col = 'b-';
    else
        col = x.color;
    end
    if  strcmpi(x.plot,'ploteye')
        ntime = -0.5:1/(GSTATE.NT-1):0.5;
        plot(ntime,Iricmat'/GSTATE.POWER(ich),col)
        xlabel('time')
        if exist('ylab','var')
            ylabel(['Eye',ylab])
        end
    elseif strcmpi(x.plot,'plotcur')
        time=0:1/GSTATE.NT:GSTATE.NSYMB-1/GSTATE.NT;
        plot(time,Iric,col);
        xlabel('time   [a.u.]');
        if exist('ylab','var')
            ylabel(['Current',ylab,'   [mA]'])
        end
    elseif strcmp(x.plot,'')
        % nothing to do
    else
        error(['the field plot must be ',...
            '''ploteye'' or ''plotcur''']);
    end
    drawnow;
end

%%%%%%%%%%%% PRINT TO FILE
if isfield(x,'print')
    if isfield(x,'ncycle')
        addnum = x.ncycle;
    else
        addnum = 0;
    end
    if isfield(x,'comp') && strcmpi(x.comp,'phase')
        x.print{1}=[x.print{1},'_phase'];
    end
    if isfield(x,'comp') && strcmpi(x.comp,'quadrature')
        x.print{1}=[x.print{1},'_quadrature'];
    end
    
    if strcmpi(x.print{2},'current')
            xax = 0:1/GSTATE.NT:GSTATE.NSYMB-1/GSTATE.NT;
            printcurrent(ich,xax,x.print,Iric,addnum);
    elseif strcmpi(x.print{2},'eye')
            xax = repmat((0:GSTATE.NT-1)/GSTATE.NT,1,GSTATE.NSYMB);
            Iric = Iric/GSTATE.POWER(ich);
            printcurrent(ich,xax,x.print,Iric,addnum);
    elseif strcmp(x.print{2},'')
            % nothing to do
    else
            error('Wrong field print');
    end    
end

%--------------------------------------------------------------------------
function printcurrent(ich,xax,name,Iric,ncycle)

%PRINTCURRENT Print the received current to file
%   PRINTCURRENT(ICH,XAX,NAME,IRIC,NCYCLE) prints the received current IRIC  
%   of channel ICH into the file NAME{1} in the global directory GSTATE.DIR. 
%   NCYCLE if ~= 0 is a counter that will be added to the output file name, 
%   and can be useful in presence of an external cycle to the function.
%   The current will be printed into a two column format where the first
%   column (the x-axis) is XAX, the second (the y-axis) is IRIC.
%
%   NAME{2} is ''eye'' or ''current'' and will be the suffix of the final
%   file. In this case IRIC is the eye pattern or the current, respectively.
%
%   The current will be printed in in a directory with the suffix ".MOD",
%
%   Author: Paolo Serena, Jan-2005
%   University of Parma, Italy

global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

MAXD = 32768;   % maximum number of points for using fprintf

Nfft = length(Iric);


if Nfft > MAXD
    types = 1;  % use save for printing (faster than fprintf)
else
    types = 0;  % use fprintf (smaller file and more readable)
end
if ~isfield(GSTATE,'DIR')
    GSTATE.DIR = tempname('.');
    mkdir(GSTATE.DIR);
    warning('optilux:eval_eye','Missing output directory in RESET_ALL. %s used',...
        GSTATE.DIR);
end

dirp = [GSTATE.DIR,'/',GSTATE.DIR,'.MOD/'];    

%%%% SET LABELS
chlab = num2str(ich,'%.3d');    %if ich < 100, adds zeros

if ncycle ~= 0
    cyclab = ['_Cycle_',num2str(ncycle,'%.4d')];    
                                % if an external cycle is present, sign it!
else
    cyclab = '';
end


%%%% WRITE TO FILE

nomewrite = [dirp,name{2},'_ch',...
        chlab,'_',name{1},cyclab,'.dat'];
printallf(nomewrite,xax,Iric,types);



