%parameters for 1D FEM
heatSource.type = 'const';                              %type of heat source field; const is const in space
heatSource.value = 0;                                   %heat source field

%boundary conditions, 1 element is left, 2 is right boundary
boundary.type = {'essential', 'natural'};
boundary.T0 = [3; 1];
boundary.q0 = [1; 1];

%Finescale conductivity params
fineCond.mu = 3;    %mean of log of lambda
fineCond.sigma = .5; %sigma of log of lambda
fineCond.nSamples = 3;

%Fine and coarse number of elements
nFine = 12;
nCoarse = 3;
assert(~mod(nFine, nCoarse), 'Error: Coarse mesh is not a divisor of fine mesh!')

%start values
theta_cf.S = 10*eye(nFine + 1);
theta_c.theta = [5; 5];
theta_c.sigma = 40;

%Define basis functions for p_c here
phi_1 = @(x) size(x, 1)/sum(1./x);
phi_2 = @(x) mean(x);
phi = {phi_1; phi_2};

%MCMC options
MCMC.method = 'MALA';                             %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 1;
MCMC.nThermalization = 5000;                              %thermalization steps
MCMC.nSamples = 50;                                    %number of samples
MCMC.nGap = 100;
MCMC.Xi_start = ones(nCoarse, 1);
%only for random walk
MCMC.MALA.stepWidth = .3;
stepWidth = 1e-1;
MCMC.randomWalk.proposalCov = stepWidth*eye(nCoarse);   %random walk proposal covariance
MCMC = repmat(MCMC, fineCond.nSamples, 1);

%stopping criterion of EM
maxIterations = 400;


%prealloc of MCMC out structure
out.samples = zeros(nCoarse, MCMC(1).nSamples);
out.log_p = zeros(MCMC(1).nSamples, 1);
out.data = cell(MCMC(1).nSamples, 1);
out.acceptance = 0;
out.log_pEnd = 0;
out = repmat(out, fineCond.nSamples, 1);

%Construct interpolation matrix W
W = zeros(nFine + 1, nCoarse);
FperC = nFine/nCoarse;
%first few lines
for i = 1:(FperC)
    W(i, 1) = (FperC - i + 1)/FperC;
    W(i, 2) = (i - 1)/FperC;
end
j = 2;
for i = 1:(nCoarse - 1)
    W((i*FperC + 1):((i + 1)*FperC), j:(j + 1)) = W(1:FperC, 1:2);
    j = j + 1;
end
W(end) = 1;





