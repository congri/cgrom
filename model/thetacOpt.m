function [F] = thetacOpt(theta_c, nData, nCoarse, sum_XNormSqMean,...
    sumPhiTXmean, sumPhiSq)
%Function returning the equation system to find the optimal theta_c

%prior derivatives; adjust prior parameters in prior functions
[~, dprior_dthetac] = log_prior_theta_c(theta_c, 'laplace');

%Equation system
F = sumPhiTXmean - sumPhiSq*theta_c + (2/(nData*nCoarse))*...
    (sum_XNormSqMean - 2*theta_c'*sumPhiTXmean + theta_c'*sumPhiSq*theta_c)*...
    dprior_dthetac;

end

