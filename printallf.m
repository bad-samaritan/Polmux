function printallf(nomewrite,xax,tempf,types)

%PRINTALLF function for fast printing to file
%   PRINTALLF(NOMEWRITE,XAX,TEMPF,TYPES) prints the vector TEMPF into the
%   file NOMEWRITE using SAVE if TYPES == 1 (faster, but bad appearance) or
%   FPRINTF if TYPES == 0 (slower).  XAX is the x-axis correspondent to
%   the y-axis TEMPF.
%
%   See also SAVE, FPRINTF
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


Nfft = length(tempf);
if types == 1
    tvar = [xax',tempf];
    eval(['save ',nomewrite,' -ASCII tvar']);
else
    fid = fopen(nomewrite,'w');   
    strpr = ['%12.5f%14.7g\n'];
    for n=1:Nfft
        fprintf(fid,strpr,xax(n),tempf(n));
    end
    fclose(fid);
end

