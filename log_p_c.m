function [log_p, d_log_p] = log_p_c(Xq, Phi, theta, sigma)
%Probabilistic mapping from fine to coarse heat conductivity
%   Xq:         Effective conductivity vector
%   x:          fine conductivities
%   phi:        basis functions
%   theta:      basis function coefficients
%   sigma:      noise
%   nFine:      Number of fine elements
%   nCoarse:    Number of coarse elements

mu  = Phi*theta;    %mean

%ignore constant prefactor
if any(Xq < 0)
    log_p = -Inf;
%     warning('Xq has negative component')
else
    log_p = - size(Xq, 1)*log(sigma) - (1/(2*sigma^2))*(Xq - mu)'*(Xq - mu);
end
if nargout > 1
    d_log_p = (1/sigma^2)*(mu - Xq);
end
    
end

