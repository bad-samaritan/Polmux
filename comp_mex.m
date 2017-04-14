function comp_mex()

%COMP_MEX compile all .c files into the directory
%   COMP_MEX is a simple routine to create .mex files starting from the .c
%   files. Please, read carefully the Matlab support guide for the options
%   of the mex function.
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


filec = dir('*.c'); % get .c files

fprintf('\n\n');
for k=1:length(filec)
    if ~exist('OCTAVE_VERSION','builtin')
        n=mex(filec(k).name);  % add here any option to the mex function
        if n == 0
            fprintf('file  %s  compiled succesfully\n',filec(k).name);
        else
            fprintf('error in compiling %s\n',filec(k).name);
        end
    else
        mex(filec(k).name);
    end
end
