function [theta_c] = optTheta_c(theta_c, fineData, nCoarse, XNormSqMean,...
    sumPhiTXmean, sumPhiSq, sumPhiSqInv, mix_sigma, prior_type, prior_hyperparam)
%% Find optimal theta_c and sigma

%% No prior on theta_c
if strcmp(prior_type, 'none')
    theta_c.theta = sumPhiSqInv*sumPhiTXmean;
    sgOpt = sigmaOpt(theta_c.theta, fineData.nSamples, nCoarse, sum(XNormSqMean),...
    sumPhiTXmean, sumPhiSq);
    theta_c.sigma = mix_sigma*theta_c.sigma + (1 - mix_sigma)*sgOpt;    
else
    %% Equation system solution options
    fsolve_options = optimoptions('fsolve');
    fsolve_options.MaxFunEvals = 3000*(numel(theta_c.theta) + 1);
    fsolve_options.MaxIter = 10000;
    fsolve_options.Algorithm = 'trust-region-dogleg';
    fsolve_options.Display = 'off';
    fsolve_options.Jacobian = 'on';
    fsolve_options.FinDiffType = 'central';     %more accurate, but slower
    %fsolve_options.DerivativeCheck = 'on';     %check derivatives with finite difference
    
    %% Generate equation system
    EqSys = @(theta) thetacOpt_Eq(theta, theta_c.theta, fineData.nSamples, nCoarse, sum(XNormSqMean),...
        sumPhiTXmean, sumPhiSq, prior_type, prior_hyperparam);
    %% Solve equation system
    [theta_c.theta, fval, exit] = fsolve(EqSys, theta_c.theta, fsolve_options)
    
    sigmaOffset = 1e-4;
    theta_c.sigma = (1 - mix_sigma)*(sigmaOpt(theta_c.theta, fineData.nSamples, nCoarse, sum(XNormSqMean),...
        sumPhiTXmean, sumPhiSq) + sigmaOffset) + mix_sigma*theta_c.sigma;
end

end

