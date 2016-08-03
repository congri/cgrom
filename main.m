%Main script for coarse-graining/reduced order modeling
clear all;

%load params
params;

Fmesh = genMesh(boundary, nFine);
Cmesh = genMesh(boundary, nCoarse);

%Generate finescale dataset
[x, Tf] = genFineData(Fmesh, heatSource, boundary, fineCond);


%get params
Cmesh.conductivity = linspace(1, 4, Cmesh.N_el);       %conductivities
params;

[Tc] = FEMmain(Cmesh, heatSource, boundary)


