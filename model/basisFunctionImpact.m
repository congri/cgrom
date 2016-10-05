function [Usq_neg, d_Usq_neg] = basisFunctionImpact(z, Xmean, new_phi, x, theta_c, PhiArray)
%% Objective function to maximize to dynamically add basis functions phi in p_c
% Input:
%   z:              params of new basis function
%   Xmean:          mean of X under qi; one column for each i
%   new_phi:        new basis function phi
%   x:              fine scale input (conductivities)
%
% Output:
%   U:              derivative of lower bound F w.r.t. theta (objective)
%   d_U:            gradient of objective function w.r.t. params z

%% Compute U
% handle to gradient of new_phi
function [d_n_phi] = d_new_phi(zin)
    [~, d_n_phi] = new_phi(zin);
end

nData = size(x, 2);
U = 0;
d_U = 0*z;

new_phi_mat = new_phi(z);
d_new_phi_mat = d_new_phi(z);
for i = 1:nData
    U = U + (Xmean(:, i)' - theta_c'*PhiArray(:, :, i)')*new_phi_mat(:, i);
    d_U = d_U + (Xmean(:, i)' - theta_c'*PhiArray(:, :, i)')*d_new_phi_mat(:, i);
end

%for maximization of modulus using fminunc
Usq_neg = -U^2;
d_Usq_neg = -2*U*d_U;

end

