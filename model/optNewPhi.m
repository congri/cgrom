function [thetaTildeOpt, zOpt] = optNewPhi(Xmean, x, theta_c, PhiArray, nFine, nCoarse)
%Find optimal parameters for new basis function phi


new_phi = @(z) phi_tilde(z, x, nFine, nCoarse);

objGrad = @(thetaTildeZ) basisFunctionGainGrad(thetaTildeZ, theta_c, Xmean, new_phi, x, PhiArray);
startValue = 4*rand(2,1) - 2;
options.Xtol = 1e-6;
options.gradTol = 1e-4;
options.initialStepSize = .001;
options.debug = true;
options.provide_objective = false;
[xOpt, gradOpt, nIter] = BFGSMaximization(objGrad, startValue, options)
thetaTildeOpt = xOpt(1);
zOpt = xOpt(2:end);









% objective = @(z) basisFunctionImpact(z, Xmean, new_phi, x, theta_c, PhiArray);
% 
% %Constrained optimization for numerical stability
% fmin_options = optimoptions('fmincon');
% fmin_options.Algorithm = 'trust-region-reflective';
% fmin_options.MaxFunctionEvaluations = 3000;
% fmin_options.MaxIter = 10000;
% fmin_options.Display = 'off';
% fmin_options.FiniteDifferenceType = 'central';     %more accurate, but slower
% fmin_options.SpecifyObjectiveGradient = true;
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% nonlcon = [];
% lb = -30;
% ub = 30;
% z0a = -2;
% z0b = 2;
% 
% [zOpta, fa] = fmincon(objective, z0a, A, b, Aeq, beq, lb, ub, nonlcon, fmin_options)
% [zOptb, fb] = fmincon(objective, z0b, A, b, Aeq, beq, lb, ub, nonlcon, fmin_options)
% if fa < fb
%     zOpt = zOpta;
% else
%     zOpt = zOptb;
% end


end

