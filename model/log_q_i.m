function [log_q, d_log_q, Tc] = log_q_i(Xi, Tf_i, theta_cf, theta_c, Phi, Fmesh, Cmesh, heatSource, boundary, W)

Cmesh.conductivity = exp(Xi);
if any(Cmesh.conductivity < 1e-6)
%     warning('very small conductivity')
%     Cmesh.conductivity
%     pause
    Cmesh.conductivity(Cmesh.conductivity < 1e-6) = 1e-6;
%     assert(~any(Cmesh.conductivity < 1e-8), 'Error: numerically 0 conductivity')
end
if any(Cmesh.conductivity > 1e6)
    Cmesh.conductivity(Cmesh.conductivity > 1e6) = 1e6;
%     warning('Very large conductivity')
end

[lg_p_c, d_lg_p_c] = log_p_c(Xi, Phi, theta_c.theta, theta_c.sigma);
[lg_p_cf, d_lg_p_cf, Tc] = log_p_cf(Tf_i, Cmesh, heatSource, boundary, W, theta_cf.S);

log_q = lg_p_cf + lg_p_c;

d_log_q = d_lg_p_c + d_lg_p_cf;
if any(imag(d_log_q))
    Xi
    theta_c.theta
    theta_c.sigma
    Phi
    d_log_q
    d_lg_p_c
    d_lg_p_cf
    error('Imaginary gradient')
end


%Finite difference gradient check
FDcheck = false;
if FDcheck
    d = 1e-4;
    gradFD = zeros(Cmesh.N_el, 1);
    for i = 1:Cmesh.N_el
        dXi = zeros(Cmesh.N_el, 1);
        dXi(i) = d;
        Cmesh.conductivity = exp(Xi) + exp(Xi).*dXi;
        
        [lg_p_c, d_lg_p_c] = log_p_c(Xi + dXi, Phi, theta_c.theta, theta_c.sigma);
        [lg_p_cf, d_lg_p_cf] = log_p_cf(Tf_i, Cmesh, heatSource, boundary, W, theta_cf.S);
        
        log_qFD = lg_p_cf + lg_p_c;
        gradFD(i) = (log_qFD - log_q)/d;
    end
    
    relgrad = gradFD./d_log_q
    
end

end

