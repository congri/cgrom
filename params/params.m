%parameters for 1D FEM
heatSource.type = 'const';                              %type of heat source field; const is const in space
heatSource.value = 0;                                   %heat source field

%boundary conditions, 1 element is left, 2 is right boundary
boundary.type = {'essential', 'natural'};
boundary.T0 = [10; 1];
boundary.q0 = [1; 400];

%Finescale conductivity params
fineCond.dist = 'uniform';  %uniform or Gaussian (dist of log conductivity)
%for Gaussian
fineCond.mu = 1.2;    %mean of log of lambda
fineCond.sigma = .3; %sigma of log of lambda
fineCond.nSamples = 150;
%for uniform
fineCond.lo = 2;
fineCond.up = 10;

%Fine and coarse number of elements
nFine = 8;
nCoarse = 4;
assert(~mod(nFine, nCoarse), 'Error: Coarse mesh is not a divisor of fine mesh!')

%start values
theta_cf.S = 1*eye(nFine + 1);
theta_cf.mu = zeros(nFine + 1, 1);
theta_c.theta = [1; 0];
theta_c.sigma = 2;

%Construct interpolation matrix W
FperC = nFine/nCoarse;
% theta_cf.W = zeros(nFine + 1, nCoarse + 1);
% %first few lines
% for i = 1:(FperC)
%     theta_cf.W(i, 1) = (FperC - i + 1)/FperC;
%     theta_cf.W(i, 2) = (i - 1)/FperC;
% end
% j = 2;
% for i = 1:(nCoarse - 1)
%     theta_cf.W((i*FperC + 1):((i + 1)*FperC), j:(j + 1)) = theta_cf.W(1:FperC, 1:2);
%     j = j + 1;
% end
% theta_cf.W(end) = 1;

%W for linear combination of coarse nodes for fine nodes in respective element
%indices for ones in matching nodes
Windices(:, 1) = (1:FperC:(nFine + 1))';
Windices(:, 2) = (1:(nCoarse + 1))';

%params indices
tempRow = 2;
tempCol = 1;
for i = 1:2*(nFine - nCoarse)
    paramIndices(i, 1) = tempRow;
    paramIndices(i, 2) = tempCol;
    if mod(i, 2)
        tempCol = tempCol + 1;
    else
        tempRow = tempRow + 1;
        tempCol = tempCol - 1;
    end
    if mod(tempRow - 1, FperC) == 0
       tempRow = tempRow + 1;
       tempCol = tempCol + 1;
    end
end
Windices = [Windices; paramIndices];
[I, J] = meshgrid(1:(nFine + 1), 1:(nCoarse + 1));
allIndices = [I(:) J(:)];
constIndices = setdiff(allIndices, paramIndices, 'rows');
clear tempRow tempCol I J allIndices;
w = rand(2*(nFine - nCoarse), 1);
%take random W as initial value
theta_cf.W = (Wmatrix(nFine, nCoarse, Windices, w) + 1) - 1;
% theta_cf.W = rand(nFine + 1, nCoarse + 1);

%Define basis functions for p_c here
%keep in mind that x is the log conductivity, x = log(lambda)
phi_1 = @(x) log(size(x, 1)/sum(1./x));
phi_2 = @(x) log(mean(x));
phi = {phi_1; phi_2};

%MCMC options
MCMC.method = 'MALA';                             %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 4;
MCMC.nThermalization = 100;                              %thermalization steps
MCMC.nSamples = 30;                                    %number of samples
MCMC.nGap = 200;
MCMC.Xi_start = zeros(nCoarse, 1);
%only for random walk
MCMC.MALA.stepWidth = 1e-3;
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
mix_S = .9;
mix_W = 0;
mix_theta = 0;





