function [phi, d_phi] = phi_tilde(z, x, nFine, nCoarse)
%new proposed phi

condTemp = exp(x);
nData = size(x, 2);
%first index is fine element number within a coarse element
%second index is coarse element number
%third index is number of data sample
for i = 1:nData
    cond(:, :, i) = reshape(condTemp(:, i), nFine/nCoarse, nCoarse);
    X(:, :, i) = reshape(x, nFine/nCoarse, nCoarse);
end

%phi is a matrix: the first index corresponds to the coarse element, the second index corresponds to
%the data point
%z.exp is a column vector
z_exp = repmat(z.exp, 1, nData);
%z.coeff is a row vector
for i = 1:nData
    phi(:, i) = - log(z.coeff*cond(:, :, i).^z_exp)';
end

end

