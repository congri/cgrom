function [Phi] = designMatrix(phi, x, nFine, nCoarse)
%Construct design matrix from basis functions phi and input data x

%every column holds fine conductivities of a single coarse element
x = reshape(x, nFine/nCoarse, nCoarse);

Phi = zeros(nCoarse, size(phi, 1));
for i = 1:nCoarse
    for j = 1:size(phi, 1)
        Phi(i, j) = phi{j}(x(:,i));
    end
end

    
    
end

