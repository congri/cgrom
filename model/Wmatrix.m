function [W, d_W] = Wmatrix(nFine, nCoarse, indices, w)
%% Computes W-matrix interpolating from coarse to fine model. Matching nodes have weight 1, fine nodes
%which are not within the corresponding coarse node have weight 0.
%       nFine:          number of fine elements
%       nCoarse:        number of coarse elements
%       w:              parameters of matrix W

assert(numel(w) == 2*(nFine - nCoarse), 'Wrong number of parameters of matrix W!')

%% Assemble W-matrix
%ones
ons = ones(nCoarse + 1, 1);
%Parameters
W = sparse(indices(:, 1), indices(:, 2), [ons; w]);

%% Compute gradient matrices
if nargout > 1
    
    d_W = zeros(nFine + 1, nCoarse + 1, 2*(nFine - nCoarse));
    for i = 1:2*(nFine - nCoarse)
        d_W(indices(nCoarse + 1 + i, 1), indices(nCoarse + 1 + i, 2), i) = 1;
    end
    
end
end

