function mdoc(func)

%MDOC Display Optilux HTML documentation in the browser
%   MDOC FUNC displays the HTML documentation for the function FUNC. 
%   
%   Under Octave the default browser is konqueror, under Matlab is the help
%   browser. Such a default behavior can be simply changed inside the
%   function.
%
%   Tip: the TAB completion under Matlab is not active for user defined
%   functions. A simple trick is to call the function as for DOC, then,
%   when the syntax is ready for a carriage return, press the Home-key 
%   (or CTRL+A) so that DOC can be changed in MDOC after digiting M.
%
%	Note: This function assumes that the Optilux documentation is in the
%		default (original) doc directory relatively to the m-files.
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


DEF_BROWSER = 'konqueror'; % default browser for Octave. Change as you like

getp = which('reset_all.m'); % get path of reset_all
pathstr = fileparts(getp);
[pathstr2, name2, ext2, versn2] = fileparts(func);

helpname = [pathstr,filesep,'..',filesep,'Documentation',...
    filesep,'optilux_doc',filesep,name2,'.html'];

if exist('OCTAVE_VERSION','builtin')
    system([DEF_BROWSER,' ',helpname,'&']);
else
    web(helpname,'-helpbrowser');
end
