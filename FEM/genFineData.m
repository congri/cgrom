function [x, Tf, PhiArray] = genFineData(domain, phi, heatSource, boundary, fineCond, nFine, nCoarse)
%Generating fine-scale dataset

x = genConductivity(fineCond, domain.N_el);

PhiArray = zeros(nCoarse, size(phi, 1), size(x, 2)); 
for i = 1:size(x, 2)
    PhiArray(:,:,i) = designMatrix(phi, x(:, i), nFine, nCoarse);
end

Tf = zeros(domain.N_el + 1, fineCond.nSamples);
for i = 1:fineCond.nSamples
    domain.conductivity = x(:,i);
    Tf(:,i) = FEMmain(domain, heatSource, boundary);
end
    
end

