function [gradient] = dF_dtheta(sigmaSqTheta_c, nData, nCoarse, sum_XNormSqMean,...
    sumPhiTXmean, sumPhiSq)
%Function returning the gradients of the lower bound F w.r.t. theta_c, sigma if there is a prior on
%theta_c and/or sigma

log_sigmaSq = sigmaSqTheta_c(1);
theta_c = sigmaSqTheta_c(2:end);

%prior derivatives; adjust prior parameters in prior functions
[~, dprior_dthetac] = log_prior_theta_c(theta_c, 'gaussian');
dprior_dsigmaSqInv = 0;

%derivative w.r.t. sigma^-2
dF_dsigmaSqInv = .5*nData*nCoarse*exp(log_sigmaSq) - .5*sum_XNormSqMean + theta_c'*sumPhiTXmean...
    - .5*theta_c'*sumPhiSq*theta_c + dprior_dsigmaSqInv;

%derivative w.r.t. theta_c
dF_dtheta_c = (1/(2*exp(log_sigmaSq)))*(sumPhiTXmean - sumPhiSq*theta_c) + dprior_dthetac;

%we need the gradient as a single vector for fsolve
gradient = [dF_dsigmaSqInv; dF_dtheta_c];

end

