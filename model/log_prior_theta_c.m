function [log_p, d_log_p] = log_prior_theta_c(theta_c, mode)
%Gives the prior log probability and derivative of parameter theta_c
% theta_c is a column vector
% mode gives the functional form of the prior

dim = size(theta_c, 1);
if strcmp(mode, 'gaussian')
    %Gaussian prior
    %hyperparameters
    mu = 0*theta_c;
    Sigma = 100*eye(dim);
    
    log_p = -.5*dim*log(2*pi) - .5*logdet(Sigma) - .5*(theta_c - mu)'*(Sigma\(theta_c - mu));
    d_log_p = - Sigma\(theta_c - mu);
elseif strcmp(mode, 'laplace')
    %Laplacian prior
    %hyperparameter
    alpha = 1e-2;
    log_p = dim*log(alpha) - dim*log(2) - alpha*sum(abs(theta_c));
    d_log_p = - alpha*sign(theta_c);
    
else
    error('Unknown prior for theta_c')
end

end

