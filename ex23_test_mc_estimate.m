% test MC_estimate
%
% In this example we study the function MC_ESTIMATE. The target is the
% measure of the mean of a vector of random variables, using for each
% element the same relative precision. This means, as we will see, that we
% cannot use the same number of samples for all elements.

clear all
clc

mc.stop = [0.1 95]; % stop criterion: 0.1 accuracy with confidence 95%
mc.method = 'mean'; % first of all we estimate the average (mean) value
mc.dim = 8;         % space dimension, i.e. estimate the mean of a vector 
                    % [1xmc.dim]
totv=100;           % run totv times the MC estimation to see if the 
                    % confidence almost coincides with the frequency
                    % interpretation of the probability.
exm = 1;            % exact average value
N=32;               % block dimension. Each run, we have N samples available

stdv = sqrt(1:mc.dim);% standard deviation of the vector. Not constant, index dependent.
totins = zeros(1,mc.dim);
for ntot=1:totv
    fprintf('ntot=%d/%d',ntot,totv);
    fprintf('\n');
    cond = ones(1,mc.dim)==1; % cond(k)==1 -> continue MC estimation for k. 
    while any(cond)
        x = exm+repmat(stdv,N,1).*randn(N,mc.dim); % Gaussian vector with std equal to stdv.
        itot = find(cond);  % search elements with low accuracy.
        for nind=itot % estimate only elements with accuracy higher than mc.stop(1)
            [cond,out] = mc_estimate(x,mc,nind);
        end
%         fprintf('%d',cond)  % Uncomment this lines to see how some
%         fprintf('\n') % conditions are met with different number of samples
    end
%     ins = stdv.^2 > out.varlim(1,:) & out.var < out.varlim(2,:);
    ins = (out.mean < exm*(1+mc.stop(1))) & (out.mean > exm*(1-mc.stop(1)));
    totins = totins+ins;
end
fprintf('\nExact average value: %.3f\n',exm);
fprintf('Example (last test):\n');
fprintf('\tEstimated average value: ');
fprintf(' %.3f',out.mean);
fprintf('\n\tNumber of runs:');
fprintf(' %d',out.nruns);
fprintf('\n\nAbout all tests:\n')
fprintf('\tNumber of success with error %.1f and confidence %.f:',...
    mc.stop(1),mc.stop(2));
fprintf('%.0f ',totins); % you will see that totins is almost mc.stop(2)
fprintf('\n')
