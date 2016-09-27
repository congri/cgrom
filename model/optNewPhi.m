function [zOpt] = optNewPhi(Xmean, x, theta_c, PhiArray, nFine, nCoarse)
%Find optimal parameters for new basis function phi


new_phi = @(z) phi_tilde(z, x, nFine, nCoarse);
objective = @(z) basisFunctionImpact(z, Xmean, new_phi, x, theta_c, PhiArray);

z0 = rand(1, nFine/nCoarse);
objective
z0
objective(z0)
[zOpt, fval] = fminunc(objective, z0)

end

