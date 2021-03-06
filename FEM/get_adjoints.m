function [lambda] = get_adjoints(K, W, Tc, y_i, mu, S, Cmesh)
%Compute adjoints for gradient computation
Wgrad = W(:, ~Cmesh.nodes); %we only need essential nodes here
% lambda = (K\Wgrad')*(S\(y_i - mu - W*Tc));

%for diagonal S
Sinv = diag(1./diag(S));
lambda = (K\Wgrad')*(Sinv*(y_i - mu - W*Tc));

end

