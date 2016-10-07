function [F, d_F] = thetacOpt_Eq(theta_c, theta_c_old, nData, nCoarse, sum_XNormSqMean,...
    sumPhiTXmean, sumPhiSq, prior_type, prior_hyperparam)
%Function returning the equation system to find the optimal theta_c

%prior derivatives; adjust prior parameters in prior functions
[~, dprior_dthetac, d2prior_d2thetac] = log_prior_theta_c(theta_c, theta_c_old, prior_type, prior_hyperparam);

%Equation system
%from sigma derivative
d_sigma = (1/(nData*nCoarse))*(sum_XNormSqMean - 2*theta_c'*sumPhiTXmean + theta_c'*sumPhiSq*theta_c);
F = (1/nData)*sumPhiTXmean - (1/nData)*sumPhiSq*theta_c + d_sigma*dprior_dthetac;

%compute the Jacobian
if nargout > 1
    d_F = - sumPhiSq;
    d_F = d_F + (1/(nData*nCoarse))*((-2*sumPhiTXmean + 2*sumPhiSq*theta_c)*dprior_dthetac');
    d_F = d_F + d_sigma*d2prior_d2thetac;
end

end

