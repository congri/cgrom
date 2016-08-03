%Main script for coarse-graining/reduced order modeling


%Generate finescale dataset
%load params
domain.N_el = 20;

params;
[x, Tf] = genFineData(domain, heatSource, boundary, fineCond)

%get params
domain.N_el = 5;
domain.conductivity = linspace(1, 4, domain.N_el);       %conductivities
params;

[Tc] = FEMmain(domain, heatSource, boundary)


