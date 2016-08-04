%Main script for coarse-graining/reduced order modeling
clear all;
addpath('./MCMCsampler')
addpath('./util')

tic;
%load params
params;

Fmesh = genMesh(boundary, nFine);
Cmesh = genMesh(boundary, nCoarse);

%Generate finescale dataset
[x, Tf] = genFineData(Fmesh, heatSource, boundary, fineCond);

theta_cf.S = eye(10);
theta_c.theta = [-1; 1];
theta_c.sigma = 1;

% If no parallel pool exists
N_Threads = 2;
if isempty(gcp('nocreate'))
    % Create with 2 workers
    parpool('local',N_Threads);
end

%store handle to every q_i in a cell array lq
lq = cell(fineCond.nSamples, 1);
parfor i = 1:fineCond.nSamples
    lq{i} = @(Xi) log_q_i(Xi, x(:,i), Tf(:,i), theta_cf, theta_c, phi, Fmesh, Cmesh, heatSource, boundary);
    %sample from every q_i
    Xi_start = [1; 1; 1];
    out(i) = MCMCsampler(lq{i}, Xi_start, MCMC);
end

runtime = toc




