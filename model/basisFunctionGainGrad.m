function [grad] = basisFunctionGainGrad(thetaTildeZ, theta_c, Xmean, phiTilde, x, PhiArray)
%Gives gradient of gain (KL_old - KL_new) of adding new basis functions in p_c
% Input:
%   z:              params of new basis function
%   Xmean:          mean of X under qi; one column for each data point i
%   new_phi:        handle to new basis function phi, output and gradient
%   x:              fine scale input (conductivities)

nData = size(x, 2);
thetaTilde = thetaTildeZ(1);
z = thetaTildeZ(2:end);
dG_dThetaTilde = 0;
dG_dz = 0*z;
[phiTilde_mat, d_phiTilde_mat] = phiTilde(z);

for i = 1:nData
    dG_dThetaTilde = dG_dThetaTilde + (Xmean(:, i)' - theta_c'*PhiArray(:, :, i)'...
        - thetaTilde*phiTilde_mat(:, i)')*phiTilde_mat(:, i);
    dG_dz = dG_dz + (Xmean(:, i)' - theta_c'*PhiArray(:, :, i)'...
        - thetaTilde*phiTilde_mat(:, i)')*thetaTilde*d_phiTilde_mat(:, i);
end

grad = [dG_dThetaTilde; dG_dz];

if(~all(isfinite(grad)))
    thetaTilde
    phiTilde_mat
    d_phiTilde_mat
    grad
    error('Non-finite gradient of G')
end
end

