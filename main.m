%Main script for coarse-graining/reduced order modeling
clear all;
addpath('./MCMCsampler')
addpath('./util')

%load params
params;

Fmesh = genMesh(boundary, nFine);
Cmesh = genMesh(boundary, nCoarse);

%Generate finescale dataset
[x, Tf] = genFineData(Fmesh, heatSource, boundary, fineCond);

Xq = [1; 1; 1];
theta = [-1; 1];
sigma = 1;
logPC = log_p_c(Xq, x(:,1), phi, theta, sigma, nFine, nCoarse)


theta_cf.S = eye(10);
theta_c.theta = [-1; 1];
theta_c.sigma = 1;
Cmesh.conductivity = [1; 1; 1];
[lq] = log_q_i(x(:,1), Tf(:,1), theta_cf, theta_c, phi, Fmesh, Cmesh, heatSource, boundary)


%get params
Cmesh.conductivity = [1; 1; 1];       %conductivities
params;

[Tc] = FEMmain(Cmesh, heatSource, boundary)

S = eye(nFine + 1);
logPCF = log_p_cf(Tf(:,1), Fmesh.coordinates, Cmesh.coordinates, Tc, S)



