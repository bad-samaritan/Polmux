function sig=chan_delay(sig,dtype,tau)

%CHAN_DELAY Add a delay to the WDM channels
%   SIG=CHAN_DELAY(SIG) adds a delay to the columns of matrix SIG generated
%   by ELECTRICSOURCE. The delay is random and uniformely distribute within
%   1 bit time.
%   The delay is positive, i.e. the signal is right-shifted on the time
%   axis.
%
%   SIG=CHAN_DELAY(SIG,DTYPE,TAU) with DTYPE='fixed' adds a deterministic
%   delay equal to TAU. In this case TAU is a vector of length equal to the
%   number of channels. In this case TAU can assume any real value.
%
%   See also ELECTRICSOURCE, RAND.
%
%   Author: Paolo Serena, 2009
%   University of Parma, Italy

global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

%Nfft = GSTATE.NSYMB*GSTATE.NT;

if (nargin == 1), tau = []; end;
Ltau = length(tau);

if nargin == 3
    if strcmp(dtype,'fixed')
        if Ltau ~= GSTATE.NCH
            error(['The delay must be a vector of length equal ',...
                'to the number of channels']);
        end
    else
        error('For a fixed delay, use the string ''fixed'' in chan_delay.m');
    end
    tauvec1 = tau;
else
    tauvec1 = rand(1,GSTATE.NCH);   
end
tauvec2 = -round(tauvec1*GSTATE.NT);    % delay in point

field_names=fieldnames(sig);
for k = 1:GSTATE.NCH
    for fn = 1:length(field_names)
        sig.(field_names{fn}) = fastshift( sig.(field_names{fn}) , -tauvec2(k) );
    end
    GSTATE.DELAY(1,k) = GSTATE.DELAY(1,k)+tauvec1(k);
end

if GSTATE.PRINT
    
    %%%% PRINT summary

    outfile = [GSTATE.DIR,'/simul_out'];

    fid = fopen(outfile,'a');

    fprintf(fid,'========================================\n');
    fprintf(fid,'===            chan_delay            ===\n');
    fprintf(fid,'========================================\n\n');

    if nargin <= 2
	fprintf(fid,'Added a random delay equal to:\n\n');
	fprintf(fid,'[');
	for k=1:GSTATE.NCH
            fprintf(fid,'%.2f ',tauvec1(k));
	end
	fprintf(fid,'] bit times x channel\n\n');
    end
    if nargin == 3
	fprintf(fid,'Added a fixed delay equal to:\n\n');
	fprintf(fid,'[');
	for k=1:GSTATE.NCH
            fprintf(fid,'%.2f ',tauvec1(k));
	end
	fprintf(fid,'] bit times x channel\n\n');
    end

    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);
       
end % end IF GSTATE.PRINT
    

