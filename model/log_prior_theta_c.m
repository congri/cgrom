function [log_p, d_log_p, d2_log_p] = log_prior_theta_c(theta_c, theta_c_old, prior_type, prior_hyperparam)
%Gives the prior log probability and derivative of parameter theta_c
% theta_c is a column vector
% mode gives the functional form of the prior

dim = size(theta_c, 1);
if strcmp(prior_type, 'gaussian')
    %Gaussian prior
    %hyperparameters
    mu = 0*theta_c;
    Sigma = prior_hyperparam;
    
    log_p = -.5*dim*log(2*pi) - .5*logdet(Sigma) - .5*(theta_c - mu)'*(Sigma\(theta_c - mu));
    d_log_p = - Sigma\(theta_c - mu);
    if nargout > 2
       d2_log_p = - inv(Sigma); 
    end
elseif strcmp(prior_type, 'laplace')
    %Laplacian prior
    %hyperparameter
    alpha = prior_hyperparam;
    log_p = dim*log(alpha) - dim*log(2) - alpha*sum(abs(theta_c));
    d_log_p = - alpha*sign(theta_c);  
    if nargout > 2
       d2_log_p = zeros(dim);
    end
elseif strcmp(prior_type, 'hierarchical')
    %Hierarchical Bayesian model with Jeffrey's hyperprior 1/Var_thetac
    %p is not a probability distribution here!
    log_p = -.5*sum((theta_c./theta_c_old).^2);
    %introduce a small offset to avoid infinities when a component of theta_c == 0
    offset = 1e-8;
    d_log_p = - theta_c./((theta_c_old + offset).^2);
    if nargout > 2
       d2_log_p = - diag(1./((theta_c_old + offset).^2)); 
    end
else
    error('Unknown prior for theta_c')
end

end

