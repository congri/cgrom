function [lq, mu] = log_q_i(Xi, Tf_i, theta_cf, theta_c, Phi, Fmesh, Cmesh, heatSource, boundary, W)
    
Cmesh.conductivity = Xi;
[Tc] = FEMmain(Cmesh, heatSource, boundary);
% mu = interp1(Cmesh.coordinates, Tc, Fmesh.coordinates)
mu = W*Tc;

lq = log_p_cf(Tf_i, mu, theta_cf.S) + log_p_c(Xi, Phi, theta_c.theta, theta_c.sigma);
    
end

