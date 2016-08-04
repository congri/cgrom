function [log_p] = log_p_c(Xq, x, phi, theta, sigma, nFine, nCoarse)
%Probabilistic mapping from fine to coarse heat conductivity
%   Xq:         Effective conductivity vector
%   x:          fine conductivities
%   phi:        basis functions
%   theta:      basis function coefficients
%   sigma:      noise
%   nFine:      Number of fine elements
%   nCoarse:    Number of coarse elements

[Phi] = designMatrix(phi, x, nFine, nCoarse);
mu  = Phi*theta;    %mean

%ignore constant prefactor
log_p = - size(Xq, 1)*log(sigma) - (1/(2*sigma^2))*(Xq - mu)'*(Xq - mu);

    
end

