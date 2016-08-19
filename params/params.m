%parameters for 1D FEM
heatSource.type = 'const';                              %type of heat source field; const is const in space
heatSource.value = 0;                                   %heat source field

%boundary conditions, 1 element is left, 2 is right boundary
boundary.type = {'essential', 'natural'};
boundary.T0 = [10; 1];
boundary.q0 = [1; 400];

%Finescale conductivity params
fineCond.mu = 3;    %mean of log of lambda
fineCond.sigma = .5; %sigma of log of lambda
fineCond.nSamples = 42;

%Fine and coarse number of elements
nFine = 6;
nCoarse = 3;
assert(~mod(nFine, nCoarse), 'Error: Coarse mesh is not a divisor of fine mesh!')

%start values
theta_cf.S = 50*eye(nFine + 1);
theta_cf.mu = zeros(nFine + 1, 1);
theta_c.theta = [.5; .5];
theta_c.sigma = 50;

%Construct interpolation matrix W
theta_cf.W = zeros(nFine + 1, nCoarse + 1);
FperC = nFine/nCoarse;
%first few lines
for i = 1:(FperC)
    theta_cf.W(i, 1) = (FperC - i + 1)/FperC;
    theta_cf.W(i, 2) = (i - 1)/FperC;
end
j = 2;
for i = 1:(nCoarse - 1)
    theta_cf.W((i*FperC + 1):((i + 1)*FperC), j:(j + 1)) = theta_cf.W(1:FperC, 1:2);
    j = j + 1;
end
theta_cf.W(end) = 1;
%take random W as initial value
%theta_cf.W = rand(nFine + 1, nCoarse + 1);

%Define basis functions for p_c here
phi_1 = @(x) size(x, 1)/sum(1./x);
phi_2 = @(x) mean(x);
phi = {phi_1; phi_2};

%MCMC options
MCMC.method = 'MALA';                             %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 88;
MCMC.nThermalization = 100;                              %thermalization steps
MCMC.nSamples = 50;                                    %number of samples
MCMC.nGap = 50;
MCMC.Xi_start = ones(nCoarse, 1);
%only for random walk
MCMC.MALA.stepWidth = .01;
stepWidth = 1e-1;
MCMC.randomWalk.proposalCov = stepWidth*eye(nCoarse);   %random walk proposal covariance
MCMC = repmat(MCMC, fineCond.nSamples, 1);

%stopping criterion of EM
maxIterations = 500;

%prealloc of MCMC out structure
out.samples = zeros(nCoarse, MCMC(1).nSamples);
out.log_p = zeros(MCMC(1).nSamples, 1);
out.data = cell(MCMC(1).nSamples, 1);
out.acceptance = 0;
out.log_pEnd = 0;
out = repmat(out, fineCond.nSamples, 1);

%Control convergence velocity - take weighted mean of adjacent parameter estimates
mix_sigma = 0;
mix_S = 0;
mix_W = 0;
mix_theta = 0;





