function [p] = p_c(Xq, x, phi, theta, sigma, nFine, nCoarse)
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

p = mvnpdf(Xq, mu, sigma*eye(length(mu)));

    
end

