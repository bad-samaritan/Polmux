function fprintfmsg(str)

%FPRINTFMSG Write a message into the file simul_out.
%   FPRINTMSG(STR) print the message STR into the file simul_out using the
%   rules of FPRINTF (\n,\t,...).
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

RIGHTLIM = 60;      % right limit for printing

if ~ischar(str)
    error('str must be a string');
end

outfile = [GSTATE.DIR,'/simul_out'];

fid = fopen(outfile,'a');
Ls = length(str);

fprintf(fid,'========================================\n');
fprintf(fid,'===          Personal Note           ===\n');
fprintf(fid,'========================================\n\n');

nmax = ceil(Ls/RIGHTLIM)-1;
if nmax == 0
    fprintf(fid,str,'\n');
else
    for k=1:ceil(Ls/RIGHTLIM)-1
        fprintf(fid,[str(1+(k-1)*RIGHTLIM:k*RIGHTLIM),'\n']);
    end
    fprintf(fid,[str(1+k*RIGHTLIM:Ls),'\n']);
end
fprintf(fid,'\n****************************************\n\n');

fclose(fid);
