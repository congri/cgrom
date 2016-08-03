function [lq] = log_q_i(x, Tf_i, theta_cf, theta_c, phi, Fmesh, Cmesh, heatSource, boundary)
    
[Tc] = FEMmain(Cmesh, heatSource, boundary);

lq = log_p_cf(Tf_i, Fmesh.coordinates, Cmesh.coordinates, Tc, theta_cf.S) + ...
    log_p_c(Cmesh.conductivity, x, phi, theta_c.theta, theta_c.sigma, Fmesh.N_el, Cmesh.N_el);

    
    
end

