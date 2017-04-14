function infoax=fieldgui(infoax)

%FIELDGUI Interactive gui for numerical propagators
%   FIELDGUI(INFOAX) Interactive propagator using matlab gui.
%   INFOAX is a structure with fields:
%
%   INFOAX.selector = operation to do (initialized to 0 before calling
%       fieldgui. After initialization, it will be controlled by fieldgui).
%   INFOAX.flag1d = [a b c]. a=1 -> plot the power of OPT. b=1 -> plot the
%       phase of OPT. c=1 -> plot the chirp of OPT. a,b,c=0 -> don't
%       plot Power,Phase,Chirp.
%   INFOAX.flag3d = same form as INFOAX.flag1d, but for 3D plot.
%   INFOAX.ch = channel number. E.g. size(OPT) = [512, 8], hence OPT
%       contains 8 channels. INFOAX.ch = [1 3] indicates that only channels 1
%       and 3 will be plotted.
%   INFOAX.dist = propagating distance.
%   INFOAX.axprop = cell array containing valid pairs of properties (see
%       axes for more details).
%       E.g. INFOAX.axprop = {'XLim',[0 10],'YLim',[0 2],'FontSize',24}
%
%   INFOAX contains also the following minor properties (created by 
%   fieldgui, not by the user):
%
%   INFOAX.numaxes = number of axes.
%   INFOAX.y3D = y-axis for 3D plot.
%   INFOAX.fig = figure handle.
%   INFOAX.editbox = vector of the button's handles.
%   INFOAX.axes = axes handles.
%
%   See Also FIBERGUI
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


persistent figint   % remember the figure between calls
% initialize GUI application
    
    gcf11plot = [0.25 0.40 0.50 0.50];  % figure dimensions with 1 plot
    gcf12plot = [0.22 0.10 0.55 0.8];    % [1,2] plot
    gcf21plot = [0.22 0.10 0.55 0.8];   % [2,1] plot
    gcf22plot = [0.10 0.10 0.80 0.80];    % [2,2] plot
    gcfhhplot = [0.002 0.07 1 0.83];     % [>2,>2] plot
    DEF_VIEW = [30 45];             % default view for 3D plot
    
    csg = length(infoax.ch);    % number of channels
    N1d = sum(infoax.flag1d);   % number of 1d plot
    N3d = sum(infoax.flag3d);   % munber of 3d plot
    numplot = csg*(N1d+N3d);    % number of plots
    switch numplot
        case 1
            gcfdim = gcf11plot; % set figure dimensions
        case 2
            gcfdim = gcf12plot;
        otherwise
            gcfdim = gcfhhplot;
    end
    ncol = ceil(numplot/2); % number of columns
    if numplot == 1
        nrow = 1;   % number of rows
    else
        nrow = ceil(numplot/8)+1;
    end
    infoax.numaxes = nrow*ncol;
    
    if  isempty(figint) 
        infoax.fig = figure('Units','Normalized','Position',gcfdim,...
            'Numbertitle','off','Toolbar','figure');  
        figint = infoax.fig;
    else
        infoax.fig = figint;
    end

    %====================================
    % information for all buttons
    
    labelColor=192/255*[1 1 1];
    top=0.95;
    bottom=0.05;
    left=0.75;
    yInitLabelPos=0.90;
    left = 0.78;
    labelWid=0.05;
    labelHt=0.06;
    btnWid = 0.08;
    btnDis = 1.2*btnWid;
    btnHt=0.07;
    % Spacing between the label and the button for the same command
    btnOffset=0.002;
    % Spacing between the button and the next command's label
    spacing=0.02;
 
    strbutton = {'Rwd','Play','Fwd','Pause'};
    
    frmBorder=0.02;

    % The boxes 
    yLabelPos=top+5.5*frmBorder;
    for nc=1:length(strbutton)
        btnNumber=nc;
           btnPos=[0.4+btnDis*(nc-3)  yLabelPos-labelHt-btnHt-btnOffset ...
               btnWid+frmBorder  btnHt];

        infoax.editbox(nc) = uicontrol( ...    % button: rwd,play,fwd,pause
        'Style','pushbutton', ...
        'Units','normalized', ...   
        'Position', btnPos, ...
        'String',strbutton(nc),...   % updsel updates infoax.selector
        'Callback',['updsel(',num2str(nc),')']);    
    end   
    btnPos=[0.4+btnDis*2  yLabelPos-labelHt-btnHt-btnOffset ...
           3*btnWid+frmBorder  btnHt];
	Minpause = 0.01;   % min value of the slider-pause
	Maxpause = 2;       % max
	Defpause = 0.5;     % default values
    
    infoax.slider = uicontrol( ...  % button slider
    'Style','slider', ...
    'Units','normalized', ...
    'Position', btnPos, ...
    'Value',Defpause,'Max',Maxpause,'Min',Minpause);
        
    %%%%%%%%%%%%%% Create axes
    ind1d3d = find([infoax.flag1d,infoax.flag3d]);
    nc = 1;
    for nch=1:csg   % channels
        for k=1:(N3d+N1d)   % axes x channel
            
        infoax.axes(nc) = subplot(nrow,ncol,nc);  % create axes
        xlabel('time [a.u.]');
        if (ind1d3d(k) == 1) ylabel(['Power ch. ',num2str(nch),' [a.u.]']); end;
        if (ind1d3d(k) == 2) ylabel(['Phase ch. ',num2str(nch),' [a.u.]']); end;
        if (ind1d3d(k) == 3) ylabel(['Chirp ch. ',num2str(nch),' [a.u.]']); end;
        if (ind1d3d(k) == 4) zlabel(['Power ch. ',num2str(nch),' [a.u.]']); end;
        if (ind1d3d(k) == 5) zlabel(['Phase ch. ',num2str(nch),' [a.u.]']); end;
        if (ind1d3d(k) == 6) zlabel(['Chirp ch. ',num2str(nch),' [a.u.]']); end;
        if (ind1d3d(k) >= 4) 
            set(infoax.axes(nc),'view',DEF_VIEW,'YLim',[1 infoax.max3d]); 
        end;
        nc = nc+1;    
        grid on;
        if isfield(infoax,'axprop') % personal settings
            if ~iscell(infoax.axprop)
                error('field axprop must be a cell array');
            end
            Laxprop = length(infoax.axprop);
            if rem(Laxprop,2)
                error('ax.prop must be of length even');
            end    
            for n=1:2:Laxprop
                set(infoax.axes(k),infoax.axprop(n),infoax.axprop(n+1));
            end
        end
        end % for k=1:(N3d+N1d)
    end     % for nch=1:csg
    selector = 2;
    set(gcf,'UserData',selector);   % selector is available between calls


