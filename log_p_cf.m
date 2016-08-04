function [log_p] = log_p_cf(Tf, Fcoord, Ccoord, Tc, S)
%Coarse-to-fine map

mu = interp1(Ccoord, Tc, Fcoord);

%ignore constant prefactor
log_p = -.5*logdet(S) - .5*(Tf - mu)'*(S\(Tf - mu));
    
    
end

