function [x, Tf, PhiArray, sumPhiSqInv] = genFineData(domain, phi, heatSource, boundary, fineCond, nFine, nCoarse)
%Generating fine-scale dataset

%% Draw conductivity/ log conductivity
if strcmp(fineCond.dist, 'uniform')
    cond = (fineCond.up - fineCond.lo)*rand(domain.N_el, fineCond.nSamples) + fineCond.lo;
    x = log(cond);
elseif strcmp(fineCond.dist, 'gaussian')
    x = normrnd(fineCond.mu, fineCond.sigma, domain.N_el, fineCond.nSamples);
    cond = exp(x);
elseif strcmp(fineCond.dist, 'binary')
    r = rand(domain.N_el, fineCond.nSamples);
    cond = fineCond.lo*ones(domain.N_el, fineCond.nSamples);
    cond(r > fineCond.p_lo) = fineCond.up;
    x = log(cond);
else
    error('unknown FOM conductivity distribution');
end

%% Compute and store design matrix for each data point
PhiArray = zeros(nCoarse, size(phi, 1), size(x, 2)); 
for i = 1:size(x, 2)
    PhiArray(:,:,i) = designMatrix(phi, cond(:, i), nFine, nCoarse);
end

%% Compute output data (finescale nodal temperatures)
Tf = zeros(domain.N_el + 1, fineCond.nSamples);
for i = 1:fineCond.nSamples
    domain.conductivity = cond(:, i);
    Tf(:,i) = FEMmain(domain, heatSource, boundary);
end

%% Compute inverse of sum_i Phi^T(x_i)^Phi(x_i)
sumPhiSq = zeros(size(phi, 1), size(phi, 1));
for i = 1:fineCond.nSamples
    sumPhiSq = sumPhiSq + PhiArray(:,:,i)'*PhiArray(:,:,i);
end
sumPhiSqInv = inv(sumPhiSq);


%% save data
save('./data/fineData/fineData', 'x', 'Tf', 'PhiArray', 'sumPhiSqInv');
    
end

