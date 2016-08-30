%Script to collect data in data arrays of EM object

%% MCMC Step width
for i = 1:fineData.nSamples
    if strcmp(MCMC(i).method, 'MALA')
        EM.MCMCStepWidth(i, k) = MCMC(i).MALA.stepWidth;
    elseif strcmp(MCMC(i).method, 'randomWalk')
        %only valid for isotropic proposals!
        EM.MCMCStepWidth(i, k) = MCMC(i).randomWalk.proposalCov(1, 1);
    elseif strcmp(MCMC(i).method, 'nonlocal')
        %do nothing
    else
        error('Unknown sampling method')
    end
end

%% Optimal params
%W matrix
EM.W(:, :, k) = theta_cf.W;
%theta
EM.theta(:, k) = theta_c.theta;
%sigma
EM.sigma(k) = theta_c.sigma;
%S
EM.S(:, :, k) = theta_cf.S;
%mu
EM.mu(:, k) = theta_cf.mu;
