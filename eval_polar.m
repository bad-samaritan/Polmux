function varargout=eval_polar(ich,x,pat,dspParameters)

%EVAL_POLAR Evaluate the polar diagram for a coherent transmission.
%   [PHASE AMPLITUDE]=EVAL_POLAR(ICH,X,PAT, DSPPARAMETERS) evaluates the 
%   polar diagram of channel ICH of a coherent transmission. PAT is the 
%   symbol pattern, after decoding (see PAT_DECODER). 
%
%   X is a structure whose fields are:
%
%   X.rec = receiver type. By now only 'coherent' is supported.
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
%   X.plot = 'y': plots the polar diagram in the active figure.
%   X.color = color string for the plot (see PLOT). E.g. 'b-'.
%   X.dpost = post compensating fiber cumulated dispersion [ps/nm]
%   X.slopez = post compensating fiber cumulated slope  [ps/nm^2]
%   X.lambda = wavelength [nm] at which the post compensating fiber
%              has a cumulated dispersion equal to X.dpost.
%   X.comp = component on which evaluate eye and calculate BER. 
%            Can be 'phase' or 'quadrature' or 'both'. In the last case the 
%            function gets two measurements over the in-phase and 
%            quadrature components, sampled with the same clock time. 
%            X.comp='both' requires PAT to be a two-column matrix with the 
%            phase/quadrature binary patterns on column 1,2, respectively.
%   X.delay = 'theory' means that the delay uses the theoretical delay 
%           saved within GSTATE.DELAY (see CREATE_FIELD). By default the 
%           delay is derived from the variance of the received samples.
%           Please note that this implementation is still very experimental
%
%   The receiver is composed of an ideal, purely linear, post compensating 
%   fiber + optical filter + coherent receiver + electrical 
%   lowpass filter. See RECEIVER_COHMIX
%   DSPPARAMETERS is a structure that describes the signal processing unit
%   of the coherent receiver. For more information see DSP4COHDEC.
%
%   PHASES and AMPLITUDES represent phase and amplitude samples of the
%   received signal.
%
%   Note 1: This function works over a copy of the electric field. All  
%           fields of the global variable GSTATE are left unchanged.
%   Note 2: All fields of X must be lowercase.
%   Note 3: Please, beware that EVAL_POLAR needs the decoded pattern, not the
%   trasmitted one. See the examples.
%
%   See also PATTERN, MYFILTER, RECEIVER_COHMIX, PAT_DECODER
%
%   Author: Marco Bertolini, 2009
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

switch x.rec
    case 'coherent'        
        [Iric,retx] = receiver_cohmix(ich,x); % get the current
        % if GSTATE.FIELDY is created along propagation, it doesn't carry
        % information, thus we must ignore it and consider only the X
        % polarization for which the pattern is defined
        if isempty(GSTATE.FIELDY_TX)
            Iric = Iric(:,1:2);
        end
    otherwise
        error('Flag X.mod must be ''coherent''');
end

if isfield(x,'delay') && strcmp(x.delay,'theory') % theoretical delay
    if isfield(x,'b2b')
        if strcmp(x.b2b,'b2b'), avgdelay = 0;end   % back-to-back: no delay
    else
        if length(GSTATE.DELAY(:,ich)) == 2
            avgdelay = 0.5*(GSTATE.DELAY(1,ich)+GSTATE.DELAY(2,ich));
        else
            avgdelay = GSTATE.DELAY(1,ich);
        end
    end
    delay = avgdelay+evaldelay(x.oftype,x.obw*0.5)+...
        evaldelay(x.eftype,x.ebw)+retx.post_delay;

else % delay estimation by cross-correlation measurement
    %     if size(Iric,2) <= 2
    %         delay(1) = corrdelay(Iric(:,1),pat(:,1),GSTATE.NT,GSTATE.NSYMB);
    %         delay(2) = corrdelay(Iric(:,2),pat(:,2),GSTATE.NT,GSTATE.NSYMB);
    %     else
    %         delay(1) = corrdelay(Iric(:,1),pat(:,1),GSTATE.NT,GSTATE.NSYMB);
    %         delay(2) = corrdelay(Iric(:,3),pat(:,3),GSTATE.NT,GSTATE.NSYMB);
    %     end
    if size(Iric,2) <= 2
        Signals = complex( Iric(:,1), Iric(:,2) );
        for a=1:GSTATE.NT
            v(a)=var(abs(Signals(a:GSTATE.NT:end)));
        end
        [val delay(1)]=min(v);
    else
        Signals = complex( Iric(:,3), Iric(:,4) );
        for a=1:GSTATE.NT
            v(a)=var(abs(Signals(a:GSTATE.NT:end)));
        end
        [val delay(2)]=min(v);
    end
    delay = delay/GSTATE.NT;
end

% apply the estimated delay
if size(Iric,2) <= 2 
    Iric = fastshift( Iric, -round(delay*GSTATE.NT) );   % Sampling in the middle of the period
else
    if size(delay,2) == 1
        delay = delay*ones(1,2);
    end
    Iric(:,1:2) = fastshift( Iric(:,1:2), -round(delay(1)*GSTATE.NT) );   % Sampling in the middle of the period
    Iric(:,3:4) = fastshift( Iric(:,3:4), -round(delay(2)*GSTATE.NT) );   % Sampling in the middle of the period
end

[phases amplitudes] = dsp4cohdec( Iric, dspParameters );

% % estimate the shift between the transmitted and the received pattern
% first_bit = abs(phases) <= pi/2;
% second_bit = phases > 0;
% 
% s1 = corrdelay(second_bit,pat(:,1),1,1024);
% s2 = corrdelay(second_bit,~pat(:,1),1,1024);
% s3 = corrdelay(second_bit,pat(:,1),1,1024);
% s4 = corrdelay(second_bit,~pat(:,1),1,1024);

varargout(1) = {phases};
varargout(2) = {amplitudes};

% EVALUATE MINIMUM HAMMING DISTANCE

% samples = amplitudes.*fastexp(phases);
% samplesr = real(samples);
% smaplesi = imag(samples);

% GRAPHICAL OUTPUT
if isfield(x,'plot') && strcmp (x.plot(1),'y')
    if ~isfield(x,'color') 
        col = 'bo';
    else
        col = x.color;
    end
    polar(phases,amplitudes,col)
    title('Modulation constellation')
    drawnow;
end
