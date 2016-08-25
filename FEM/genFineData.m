function [x, Tf, PhiArray] = genFineData(domain, phi, heatSource, boundary, fineCond, nFine, nCoarse)
%Generating fine-scale dataset

%x = genConductivity(fineCond, domain.N_el);
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

PhiArray = zeros(nCoarse, size(phi, 1), size(x, 2)); 
for i = 1:size(x, 2)
    PhiArray(:,:,i) = designMatrix(phi, cond(:, i), nFine, nCoarse);
end

Tf = zeros(domain.N_el + 1, fineCond.nSamples);
for i = 1:fineCond.nSamples
    domain.conductivity = cond(:, i);
    Tf(:,i) = FEMmain(domain, heatSource, boundary);
end

save('./data/fineData/fineData', 'x', 'Tf', 'PhiArray');
    
end

