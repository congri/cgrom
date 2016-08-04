function [log_p] = log_p_cf(Tf, mu, S)
%Coarse-to-fine map

%ignore constant prefactor
% log_p = -.5*logdet(S, 'chol') - .5*(Tf - mu)'*(S\(Tf - mu));
%diagonal S
log_p = -.5*sum(log(diag(S))) - .5*(Tf - mu)'*(S\(Tf - mu));
    
end

