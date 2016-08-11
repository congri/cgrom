function [log_p, d_log_p, Tc] = log_p_cf(Tf, Cmesh, heatSource, boundary, W, S)
%Coarse-to-fine map
%ignore constant prefactor
%log_p = -.5*logdet(S, 'chol') - .5*(Tf - mu)'*(S\(Tf - mu));
%diagonal S

[Tc, d_r, K] = FEMmain(Cmesh, heatSource, boundary);
WTc = W*Tc;

%only for diagonal S!
assert(isdiag(S), 'Error: matrix S not diagonal');
Sinv = diag(1./diag(S));
log_p = -.5*sum(log(diag(S))) - .5*(Tf - WTc)'*(Sinv*(Tf - WTc));

if nargout > 1
    lambda = get_adjoints(K, W, Tc, Tf, 0*Tf, S, Cmesh);
    d_log_p = - d_r'*lambda;
    
    %Finite difference gradient check
    FDcheck = false;
    d = 1e-6;
    if FDcheck
        d_log_pFD = zeros(Cmesh.N_el, 1);
        CmeshFD = Cmesh;
        for i = 1:Cmesh.N_el
            dXq = zeros(Cmesh.N_el, 1);
            dXq(i) = d;
            CmeshFD.conductivity = Cmesh.conductivity + dXq;
            TcFD = FEMmain(CmeshFD, heatSource, boundary);
            muFD = W*TcFD;
            log_pFD = -.5*sum(log(diag(S))) - .5*(Tf - muFD)'*(Sinv*(Tf - muFD));
            d_log_pFD(i) = (log_pFD - log_p)/d;
        end 
        
        relGrad = d_log_pFD./d_log_p
        d_log_p
        d_log_pFD

        if any(relGrad > 1.2) || any(relGrad < .8)
            relGrad
            pause
        end
        
    end
end

end

