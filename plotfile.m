function plotfile(file)

%PLOTFILE Plot file from disk.
%   PLOTFILE(FILE) plot the double-column data stored in the file
%   named FILE from disk. FILE can be a compressed file (with 
%   extension .gz,.bz2,.zip).
%
%   See also LOADFILE
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


file_t = loadfile(file);
[a,b] = size(file_t);
if b == 2
    plot(file_t(:,1),file_t(:,2))
else
    error('FILE must be a double-column file');
end

% --------------------------------------------------------------------


function y=loadfile(file)

%LOADFILE Load file (possibly compressed) from disk.
%   Y=LOADFILE(FILE) load the file of data named FILE. It coincide
%   with the function LOAD except that FILE can be a compressed
%   file (with extension .gz,.bz2,.zip), provided that
%   the correspondent uncompressor tool is installed in your system.
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


[pathstr,namestr,lastpart,nameversn] = fileparts(file);

exO = exist('OCTAVE_VERSION');
if exO
    name_t = tmpnam('','tp');    % generate a unique temp name
else
    name_t = tempname;
end
switch lastpart
    case '.gz'
        system(['gunzip ',file,' -c > ',name_t]);
        y = load(name_t);
        if exO
            system(['rm ',name_t]);
        else
            delete(name_t);
        end
    case '.bz2'
        system(['bunzip2 ',file,' -c > ',name_t]);
        y = load(name_t);        
        if exO
            system(['rm ',name_t]);
        else
            delete(name_t);
        end            
    case '.zip'
        system(['unzip -p ',file,' > ',name_t]);
        y = load(name_t);        
        if exO
            system(['rm ',name_t]);
        else
            delete(name_t);
        end            
    otherwise
        y = load(file);        
end
