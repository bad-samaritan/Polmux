function [best_eo,best_dpost] = best_eye(ich,x,pat)

%BEST_EYE Search algorithm for the best eye opening.
%   BEYE=BEST_EYE(ICH,X,PAT) searches the minimum eye closure penalty BEYE 
%   of a  non-coherent transmission by varying the post compensating fiber 
%   dispersion in front of the receiver (E.g. see RECEIVER_OOK or 
%   RECEIVER_DPSK). 
%
%   [BEYE,BPOST]=BEST_EYE(ICH,X,PAT) also returns in BPOST the best post 
%   compensating fiber cumulated dispersion in [ps/nm]. 
%
%   The eye closure penalty BEYE is defined as:
%
%   BEYE = 10*log10(eob2b/eo1)   [dB]
%
%   where eo1=max(min1-max0) is the eye opening after propagation, being
%   min1 and max0 the worst mark/space samples. eob2b is the eye opening in
%   back to back configuration. 
%
%   PAT contains the symbols pattern after decoding (see PAT_DECODER).
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
%   X.dpost =  post compensating fiber cumulated dispersion [ps/nm]: Vector
%              of two elements [d1 d2] that is the range within it the best 
%              post compensation is searched. X.dpost can also be a scalar:
%              in this case, the function returns the eye closure penalty 
%              at this value of X.dpost.
%   X.slopez = post compensating fiber cumulated slope  [ps/nm^2]
%   X.lambda = wavelength [nm] at which the post compensating fiber
%              has a cumulated dispersion equal to X.dpost.
%
%   Optional parameters for X:
%
%   X.tol =    tolerance of the golden search algorithm (fractional 
%              precision: +/-tol).
%   X.print=   structure for print. E.g. X.print = {'nomefile','eye'} or
%              X.print = {'nomefile','current'}, prints to file nomefile  
%              the eye or the current, respectively. nomefile will be
%              placed into GSTATE.DIR within a directory ending with '.MOD'.
%   X.plot =   'ploteye': plots the eye in the active figure; 'plotcur' 
%              plots the received current.
%   X.color=   color string for the plot (see plot.m). E.g. 'b-'.
%   X.ts =     Fixed sampling time (-0.5 <= X.ts <= 0.5). 
%   X.comp = component on which evaluate eye and calculate BER (dqpsk 
%            modulation only). Can be 'phase' or 'quadrature' or 'both'. In
%            the last case the function gets two measurements over the
%            in-phase and quadrature components, sampled with the same
%            clock time. X.comp='both' requires PAT to be a two-column
%            matrix with the phase/quadrature binary patterns on column 1,2,
%            respectively.
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
%   The best eye opening is searched through a golden search algorithm (see
%   Numerical Recipes in Fortran 77, W. H. Press, S. A. Teukolsky and 
%   W. T. Vetterling, University of Cambridge.)
%
%   Note 1: The golden search algorithm works when only one minimum
%         is present within the range [d1 d2]. Otherwise the returned 
%         BEYE is just one of the local min, and may not be the lowest.
%   Note 2: This function works over a copy of the electric field. All  
%         fields of the global variable GSTATE are left unchanged.
%
%   See also EVAL_EYE, BEST_SP, PATTERN, MYFILTER, RECEIVER_OOK,  
%            RECEIVER_DPSK, RECEIVER_DQPSK, CORRDELAY
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

if ~isfield(x,'tol') 
    tol = 1e-6;    % default tolerance
else
    tol = x.tol;
end; 
Ldp = length(x.dpost);
if Ldp == 1
    best_dpost = x.dpost;
    [best_eo,ts]=eval_eye(ich,z,pat);
    z.plot = '';
    z.print = {'',''};
elseif Ldp == 2
    if isfield(z,'plot') || isfield(z,'print')
        warning('optilux:best_eye',...
            'Plot or print turned off during the search of the optimal post');
        z.plot = '';
        z.print = {'',''};
    end
    [best_eo,best_dpost,x12] = goldeneye(ich,pat,z,tol);
else
    error('The field dpost must be of 1 (fixed post) or 2 elements (range of search)');
end



%%%%%%%%%%%%%%%%%%%% eye opening in back-to-back
z.b2b = 'b2b';
z.dpost = 0;
[eob2b,tsb2b] = eval_eye(ich,z,pat);

best_eo = 10*log10(eob2b/best_eo);

if (nargout >= 1), varargout(1)={best_eo}; end
if (nargout == 2), varargout(2)={best_dpost}; end

if GSTATE.PRINT

    %%%% PRINT summary

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===            best_eye              ===\n');
    fprintf(fid,'========================================\n\n');
    
    if isfield(x,'comp') 
        fprintf(fid,'%s Receiver (%s):\n\n',upper(x.rec),x.comp);
    else
        fprintf(fid,'%s Receiver:\n\n',upper(x.rec));
    end
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
	fprintf(fid,'Eye closure penalty: %.4f [dB]\n',best_eo);
    else
	fprintf(fid,'Best eye closure penalty %.4f [dB]\n',best_eo);
    end

    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT

%--------------------------------------------------------------------------
function [best_eo,best_dpost,x12] = goldeneye(ich,pat,z,tol)

% Golden search algorithm for the optimal post compensation.

R=0.61803399;
C=1-R;

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
[f1,ts]=eval_eye(ich,z,pat);    % The initial function evaluations. 
                         % Note that we never need to
                         % evaluate the function at the original endpoints
z.dpost = x2;                         
[f2,ts]=eval_eye(ich,z,pat);
while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) % while loop: we keep returning here.
    if(f2>f1)    %One possible outcome,
        x0=x1;   %its housekeeping,
        x1=x2;
        x2=R*x1+C*x3;
        f1=f2;
        z.dpost = x2;
        [f2,ts]=eval_eye(ich,z,pat);   % and a new function evaluation.
    else            %The other outcome,
        x3=x2;
        x2=x1;
        x1=R*x2+C*x0;
        f2=f1;
        z.dpost = x1;
        [f1,ts]=eval_eye(ich,z,pat);     % and its new function evaluation.
    end
end
if(f1>f2) %We are done. Output the best of the two current values.
    best_eo=f1;
    best_dpost=x1;
else
    best_eo=f2;
    best_dpost=x2;
end
x12 = abs(x1)+abs(x2);
