function [x, Tf] = genFineData(domain, heatSource, boundary, fineCond)
%Generating fine-scale dataset

x = genConductivity(fineCond, domain.N_el);

Tf = zeros(domain.N_el + 1, fineCond.nSamples);
for i = 1:fineCond.nSamples
    domain.conductivity = x(:,i);
    Tf(:,i) = FEMmain(domain, heatSource, boundary);
end
    
end

