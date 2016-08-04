function [lq] = log_q_i(Xi, x, Tf_i, theta_cf, theta_c, phi, Fmesh, Cmesh, heatSource, boundary)
    
Cmesh.conductivity = Xi;
[Tc] = FEMmain(Cmesh, heatSource, boundary);

lq = log_p_cf(Tf_i, Fmesh.coordinates, Cmesh.coordinates, Tc, theta_cf.S) + ...
    log_p_c(Xi, x, phi, theta_c.theta, theta_c.sigma, Fmesh.N_el, Cmesh.N_el);

    
    
end

