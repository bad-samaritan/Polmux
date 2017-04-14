function varargout = best_sp(ich,x,pat)

%BEST_SP Search algorithm for the best OSNR penalty vs. back-to-back.
%   VARARGOUT=BEST_SP(ICH,X,PAT) searches the minimum osnr penalty, 
%   also called sensitivity penalty, of channel ICH  vs. back to 
%   back transmission, by varying the post compensating fiber dispersion. 
%
%   The function works for non-coherent transmissions. The target bit error 
%   rate (ber) of the osnr penalty is measured through the Karhunen-Loeve 
%   (kl) method (see BER_KL for more details). 
%
%   VARARGOUT can assume a variable number of arguments, from 1 to 4. In 
%   the complete case it is 
%
%       VARARGOUT=[BSP,BPOST,WARN,ECP] 
%
%   where BSP is the best OSNR penalty [dB] using the best post-dispersion 
%   BPOST [ps/nm]. WARN is a flag equal to 1 when the function found a 
%   sensitivity penalty equal to NaN during the search of the optimal post, 
%   that corresponds to a possible wrong result. ECP is the eye closure 
%   penalty [dB] using BPOST.
%
%   PAT is the symbols pattern. X is a structure whose fields are:
%
%   X.rec = receiver type. Valid arguments are: 'ook','psbt','nf-dpsk' to 
%           use RECEIVER_OOK, 'dpsk' to use RECEIVER_DPSK, 'dqpsk' to use 
%           with RECEIVER_DQPSK.
%   X.ber =    target ber at which the algorithm measures the OSNR penalty.
%   X.oftype = optical filter type (see MYFILTER)
%   X.obw = optical filter 3 dB bandwidth normalized to the symbol rate
%   X.oord = optical filter order (if X.oftype is 'supergaussian')
%   X.eftype = electrical filter type (see MYFILTER)
%   X.ebw = electrical filter 3-dB bandwidth normalized to the symbol rate
%   X.eord = electrical filter order (if X.eftype is 'supergaussian')
%   X.dpost =  post compensating fiber cumulated dispersion [ps/nm]: Vector
%              of two elements [d1 d2] that is the range within it the best 
%              post compensation is searched. X.dpost can also be a scalar:
%              in this case, the function returns the eye closure penalty 
%              at this value of X.dpost.
%   X.slopez = post compensating fiber cumulated slope  [ps/nm^2]
%   X.lambda = wavelength [nm] at which the post compensating fiber
%              has a cumulated dispersion equal to X.dpost.
%   X.osnr =   Optical signal to noise ratios (osnr), [dB], over which the 
%              ber is evaluated. The osnr is over a conventional bandwidth 
%              of 0.1 nm and is measured immediately before the receiver. 
%              X.osnr refers to X.poln noise polarizations. X.osnr is a
%              vector.
%   X.poln =   Noise polarizations, 1 or 2. Note: X.poln is independent
%              from the signal polarizations, e.g. the algorithm can work 
%              with two noise polarizations and just one signal 
%              polarization.
%
%   X have also the following parameters required by the kl-method:
%
%   X.eta =    Bandwidth expansion factor. The kl method samples the signal
%              and the noise up to a frequency equal to X.eta times the 
%              bandwidth of the optical filter. 
%              Usually is (faster) 1 < X.eta < 3 (slower).
%   X.mu =     Time expansion factor. The memory of the receiver is X.mu 
%              times the time duration of the memory devices inside the 
%              receiver, i.e. the optical/electrical filter. For DPSK there 
%              is an additional memory due to the Mach-Zehnder delay 
%              interferometer. The memory of such devices is approximated 
%              by the inverse of their bandwidth, as suggested in [1].
%              Usually is (faster) 1 < X.mu < 10 (slower).
%   X.saddle = 'yes': the ber is evaluated through the saddle point
%              approximation (faster). 'no': the ber is evaluated by 
%              numerical integration of the moment generating function 
%              (slower, but more accurate).
%
%   For more details about X.eta, X.mu and X.saddle see BER_KL.
%   X can also have the optional parameters:
%
%   X.interp = interpolation method for finding X.ber, see INTERP1.
%              Default is 'spline'.
%   X.extrap = 'yes': X.ber can be extrapolated outside X.osnr, see
%              INTERP1. 'no': If X.ber is outside the range X.osnr the 
%              function returns VARARGOUT(2) = NaN, which is also the 
%              default strategy.
%   X.tol =    tolerance of the golden search algorithm (fractional 
%              precision: +/-tol).
%   X.plot =   'ploteye': plots the eye in the active figure; 'plotcur' 
%              plots the received current.
%   X.color =  color string for the plot (see PLOT). E.g. 'b-'.
%   X.print =  structure for print. E.g. X.print = {'nomefile','eye'} or
%              X.print = {'nomefile','current'}, prints to file nomefile  
%              the eye or the current, respectively. nomefile will be place 
%              into GSTATE.DIR within a directory ending with '.MOD'.
%   X.delay = 'theory' means that the delay uses the theoretical delay 
%           saved within GSTATE.DELAY (see CREATE_FIELD). By default the 
%           delay is measured by a cross-correlation measurement between 
%           the received current and an artificial pulse amplitude 
%           modulation (PAM) signal with ideal non-return to zero bits with
%           symbols equal to PAT. The correlation method is useful in 
%           presence of polarization mode dispersion (PMD).
%
%   X.ts =     Fixed sampling time (-0.5 <= X.ts <= 0.5). 
%
%   The receiver is composed of an ideal, purely linear, post compensating 
%   fiber + optical filter + optical to electrical converter + electrical 
%   lowpass filter. For example see RECEIVER_OOK or RECEIVER_DPSK.
% 
%   The best ber is searched through a golden search algorithm (see
%   Numerical Recipes in Fortran 77, W. H. Press, S. A. Teukolsky and 
%   W. T. Vetterling, University of Cambridge.)
%
%   Note 1: The golden search algorithm works when only one minimum
%         is present within the range [d1 d2]. Otherwise the returned 
%         BSP is just one of the local minima, and may not be the lowest.
%   Note 2: This function works over a copy of the electric field. All  
%         fields of the global variable GSTATE are left unchanged.
%
%   See also PATTERN, BER_KL, BEST_EYE, MYFILTER, RECEIVER_OOK, 
%            RECEIVER_DPSK, RECEIVER_DQPSK, EVAL_EYE, CORRDELAY
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


z = x;  % copy
memprint = GSTATE.PRINT;
GSTATE.PRINT=false; % turn off printing during the search

if ~isfield(x,'tol') 
    tol = 1e-6;    % default tolerance
else
    tol = x.tol;
end; 
Ldp = length(x.dpost);
if Ldp == 1
    best_dpost = x.dpost;
    [ts,best_osnr,best_eo]=ber_kl(ich,z,pat);
    warn = 0;
    z.plot = '';
    z.print = {'',''};
elseif Ldp == 2
    if isfield(z,'plot') || isfield(z,'print')
        warning('optilux:best_sp',...
            'Plot or print turned off during the search of the optimal post');
        z.plot = '';
        z.print = {'',''};
    end
    [best_osnr,best_dpost,warn,x12,best_eo] = goldenber(ich,pat,z,tol);
else
    error('the field dpost must be of 1 or 2 elements');
end

%%%%%%%%%%%%%%%%%%%% OSNR in back-to-back
z.b2b = 'b2b';
z.dpost = 0;
[ts,osnr_b2b,eob2b] = ber_kl(ich,z,pat);

best_sp = best_osnr-osnr_b2b;   % OSNR or sensitivity penalty
best_ecp = 10*log10(eob2b/best_eo);  % Eye closure penalty [dB]

if warn
    warning('optilux:best_sp','Found NaN during the search of the best post');
end
if (nargout >= 1), varargout(1)={best_sp}; end
if (nargout >= 2), varargout(2)={best_dpost}; end
if (nargout >= 3), varargout(3)={warn}; end
if (nargout == 4), varargout(4)={best_ecp}; end

GSTATE.PRINT=memprint; % recover the initial state

if GSTATE.PRINT

    %%%% PRINT summary

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===            best_sp               ===\n');
    fprintf(fid,'========================================\n\n');

    fprintf(fid,'%s Receiver:\n\n',upper(x.rec));
    if Ldp == 1
	fprintf(fid,'Post compensating cumulated dispersion:\n');
	fprintf(fid,['  %.4f [ps/nm]@ lambda=%.2f',...
            '[nm] (cum. slope = %.2g [ps/nm^2])\n'],x.dpost,x.lambda,x.slopez);
    else
	fprintf(fid,'Post compensating fiber inital range:\n');
	fprintf(fid,'[%.3f %.3f] [ps/nm]\n',x.dpost(1),x.dpost(2));
	fprintf(fid,'post values @lambda = %.2f (cum slope = %.2g [ps/nm^2])\n\n',...
            x.lambda,x.slopez);
	fprintf(fid,'Best post: %.4f [ps/nm] (+/-: %.2g [ps/nm])\n\n',...
            best_dpost,tol*x12);
    end
    fprintf(fid,'Optical filter: %s of bandwidth %.3f\n',x.oftype,x.obw);
    fprintf(fid,'Electrical filter: %s of bandwidth %.3f\n',x.eftype,x.ebw);
    if Ldp == 1
	fprintf(fid,'OSNR penalty: %.4f [dB] (OSNR b2b = %.4f [dB])\n',...
            best_sp,osnr_b2b);

    else
	fprintf(fid,'Best OSNR penalty %.4f [dB] (OSNR b2b = %.4f [dB])\n',...
            best_sp,osnr_b2b);
    end
    fprintf(fid,'Eye Closure penalty: %.4f [dB]\n',best_ecp);
    if warn
	fprintf(fid,'\nFound NaN during the search\n');
    end

    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT

%--------------------------------------------------------------------------
function [best_osnr,best_dpost,warn,x12,best_eo] = goldenber(ich,pat,z,tol)

% Golden search algorithm for the optimal post compensation.

R=0.61803399;
C=1-R;
warn = 0;

ax = z.dpost(1);        % initialization: left corner
bx = 0.5*(z.dpost(1)+z.dpost(2));   % first guess
cx = z.dpost(2);                    % right corner

x0=ax;  % start golden search
x3=cx;
if(abs(cx-bx)>abs(bx-ax)) % Make x0 to x1 the smaller segment,
    x1=bx;
    x2=bx+C*(cx-bx);                 % and �ll in the new point to be tried.
else
    x2=bx;
    x1=bx-C*(bx-ax);
end
z.dpost = x1;
[ts,f1,e1]=ber_kl(ich,z,pat);    % The initial function evaluations. 
                         % Note that we never need to
                         % evaluate the function at the original endpoints
if isnan(f1), warn = 1; end;                         
z.dpost = x2;                         
[ts,f2,e2]=ber_kl(ich,z,pat);
if isnan(f2), warn = 1; end;                         

while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) % while loop: we keep returning here.
    if(f2<f1)    %One possible outcome,
        x0=x1;   %its housekeeping,
        x1=x2;
        x2=R*x1+C*x3;
        f1=f2;
        z.dpost = x2;
        [ts,f2,e2]=ber_kl(ich,z,pat);   % and a new function evaluation.
        if isnan(f2), warn = 1; end;                         
    else            %The other outcome,
        x3=x2;
        x2=x1;
        x1=R*x2+C*x0;
        f2=f1;
        z.dpost = x1;
        [ts,f1,e1]=ber_kl(ich,z,pat);     % and its new function evaluation.
        if isnan(f1), warn = 1; end;                         
    end
end
if(f1<f2) %We are done. Output the best of the two current values.
    best_osnr=f1;
    best_dpost=x1;
    best_eo=e1;
else
    best_osnr=f2;
    best_dpost=x2;
    best_eo=e2;
end
x12 = abs(x1)+abs(x2);
