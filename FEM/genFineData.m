function [x, Tf] = genFineData(domain, heatSource, boundary, fineCond)
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

%% Compute output data (finescale nodal temperatures)
Tf = zeros(domain.N_el + 1, fineCond.nSamples);
for i = 1:fineCond.nSamples
    domain.conductivity = cond(:, i);
    Tf(:,i) = FEMmain(domain, heatSource, boundary);
end
    
end

