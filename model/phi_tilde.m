function [phi, d_phi] = phi_tilde(z, x, nFine, nCoarse)
%new proposed phi

condTemp = exp(x);
nData = size(x, 2);
finePerCoarse = nFine/nCoarse;
%first index is fine element number within a coarse element
%second index is coarse element number
%third index is number of data sample
for i = 1:nData
    cond(:, :, i) = reshape(condTemp(:, i), finePerCoarse, nCoarse);
    X(:, :, i) = reshape(x(:, i), finePerCoarse, nCoarse);
end

%phi is a matrix: the first index corresponds to the coarse element, the second index corresponds to
%the data point
if(z == 0)
%     warning('z == 0, what is the value of the derivative?')
    sum_exp = sum(exp(X), 1);
    for k = 1:nCoarse
        for i = 1:nData
            phi(k, i) = (1/finePerCoarse)*sum_exp(1, k, i);
            d_phi(k, i) = 0;    %this is wrong
        end
    end
elseif z < -30
%     warning('z < -100, what is the value of the derivative?')
    for k = 1:nCoarse
        for i = 1:nData
            phi(k, i) = log(min(X(:, k, i)));
            d_phi(k, i) = 0;    %this is wrong
        end
    end
elseif z > 30
%     warning('z > 100, what is the value of the derivative?')
    for k = 1:nCoarse
        for i = 1:nData
            phi(k, i) = log(max(X(:, k, i)));
            d_phi(k, i) = 0;    %this is wrong
        end
    end
else
    sum_exp = sum(exp(z*X), 1);
    sum_exp_x = sum(X.*exp(z*X), 1);
    for k = 1:nCoarse
        for i = 1:nData
            phi(k, i) = (1/z)*(-log(finePerCoarse) + log(sum_exp(1, k, i)));
            d_phi(k, i) = (1/z)*(sum_exp_x(1, k, i)/sum_exp(1, k, i) - phi(k, i));
        end
    end
end

end

