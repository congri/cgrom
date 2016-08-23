function [log_p, d_log_p, Tc] = log_p_cf(Tf_i, Cmesh, heatSource, boundary, W, S)
%Coarse-to-fine map
%ignore constant prefactor
%log_p = -.5*logdet(S, 'chol') - .5*(Tf - mu)'*(S\(Tf - mu));
%diagonal S

[Tc, d_r, K] = FEMmain(Cmesh, heatSource, boundary);
%d_r is derivative w.r.t. conductivity lambda. we want derivative w.r.t. x = log(lambda).
d_rx = d_r*diag(Cmesh.conductivity);

WTc = W*Tc;

%only for diagonal S!
assert(isdiag(S), 'Error: matrix S not diagonal');
Sinv = diag(1./diag(S));
log_p = -.5*sum(log(diag(S))) - .5*(Tf_i - WTc)'*(Sinv*(Tf_i - WTc));

if nargout > 1
    lambda = get_adjoints(K, W, Tc, Tf_i, 0*Tf_i, S, Cmesh);
    d_log_p = - d_rx'*lambda;
    
    %Finite difference gradient check
    FDcheck = false;
    d = 1e-3;
    if FDcheck
        d_log_pFD = zeros(Cmesh.N_el, 1);
        CmeshFD = Cmesh;
        for i = 1:Cmesh.N_el
            dXq = zeros(Cmesh.N_el, 1);
            dXq(i) = d;
            CmeshFD.conductivity = Cmesh.conductivity + Cmesh.conductivity.*dXq;
            TcFD = FEMmain(CmeshFD, heatSource, boundary);
            muFD = W*TcFD;
            log_pFD = -.5*sum(log(diag(S))) - .5*(Tf_i - muFD)'*(Sinv*(Tf_i - muFD));
            d_log_pFD(i) = (log_pFD - log_p)/d;
        end 
        
        relGrad = d_log_pFD./d_log_p
%         d_log_p
%         d_log_pFD

        if any(relGrad > 1.2) || any(relGrad < .8)
            relGrad
            pause
        end
        
    end
end

end

