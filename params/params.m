%% parameters for 1D FEM
heatSource.type = 'const';                              %type of heat source field; const is const in space
heatSource.value = 0;                                   %heat source field

%% boundary conditions, 1. element is left, 2. is right boundary
boundary.type = {'essential', 'natural'};
boundary.T0 = [10; 1];
boundary.q0 = [1; 400];

%% Finescale conductivity params
fineData.genData = true;
fineData.dist = 'binary';  %uniform, gaussian or binary (dist of log conductivity)
fineData.nSamples = 16;
if strcmp(fineData.dist, 'gaussian')
    fineData.mu = 1.2;    %mean of log of lambda
    fineData.sigma = .3; %sigma of log of lambda
elseif (strcmp(fineData.dist, 'uniform') || strcmp(fineData.dist, 'binary'))
    %for uniform & binary
    fineData.lo = 2;
    fineData.up = 10;
    contrast = fineData.up/fineData.lo;
    %for binary
    if strcmp(fineData.dist, 'binary')
        fineData.p_lo = .4;
    end
else
    error('unknown fineCond distribution');
end

%% Fine and coarse number of elements
nFine = 32;
nCoarse = 4;
assert(~mod(nFine, nCoarse), 'Error: Coarse mesh is not a divisor of fine mesh!')

%% Coarse to fine interpolation matrix W
FperC = nFine/nCoarse;
Winterp = true;
if Winterp
    theta_cf.W = zeros(nFine + 1, nCoarse + 1);
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
else
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
end
% theta_cf.W = rand(nFine + 1, nCoarse + 1);

%% Basis functions for p_c
%keep in mind that x is the log conductivity, x = log(lambda)
phi_1 = @(x) log(size(x, 1)/sum(1./x));
phi_2 = @(x) log(mean(x));
phi_3 = @(x) log(max(x));
phi_4 = @(x) log(min(x));
phi_5 = @(x) log(mean(x.^2));
phi_6 = @(x) log(mean(x.^3));
phi_7 = @(x) log(mean(x.^4));
phi_8 = @(x) log(mean(x.^5));
phi_9 = @(x) log(secOrderFeature(x));   %sum_i x_i*x_{i + 1}
phi_10 = @(x) log(x(1));
phi_11 = @(x) log(x(2));
phi_12 = @(x) log(x(3));
phi_13 = @(x) log(x(4));
phi_14 = @(x) log(x(1))*log(x(4));
phi_15 = @(x) mean(log(x));             %geometric mean

phi = {phi_15};
nBasis = numel(phi);

%% start values
theta_cf.S = 1*eye(nFine + 1);
theta_cf.mu = zeros(nFine + 1, 1);
theta_c.theta = (1/size(phi, 1))*ones(size(phi, 1), 1);
theta_c.sigma = .5;
%what kind of prior for theta_c
prior_type = 'hierarchical';    %hierarchical, laplace, gaussian or none
%prior hyperparams; obsolete for no prior and hierarchical prior
% prior_hyperparam = 100*eye(size(phi, 1));     %variance of prior gaussian
prior_hyperparam = .0001;                       %Exponential decay parameter for laplace



%% MCMC options
MCMC.method = 'MALA';                             %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 13;
%thermalization steps; we perform a deterministic optimization of qi before each sampling, hence
%there is no thermalization needed
MCMC.nThermalization = 0;
MCMC.nSamples = 50;                                    %number of samples
MCMC.nGap = 100;
MCMC.Xi_start = mvnrnd(zeros(nCoarse, 1), theta_c.sigma);
%only for random walk
MCMC.MALA.stepWidth = 3e-5;
stepWidth = 1e-3;
MCMC.randomWalk.proposalCov = stepWidth*eye(nCoarse);   %random walk proposal covariance
MCMC = repmat(MCMC, fineData.nSamples, 1);

%% MCMC options to find step width
MCMCstepWidth = MCMC;
for i = 1:fineData.nSamples
    MCMCstepWidth(i).nSamples = 2;
    MCMCstepWidth(i).nGap = 300;
end


%prealloc of MCMC out structure
out.samples = zeros(nCoarse, MCMC(1).nSamples);
out.log_p = zeros(MCMC(1).nSamples, 1);
out.data = cell(MCMC(1).nSamples, 1);
out.acceptance = 0;
out.log_pEnd = 0;
out = repmat(out, fineData.nSamples, 1);

%% EM options
%Control convergence velocity - take weighted mean of adjacent parameter estimates
mix_sigma = 0;
mix_S = 0;
mix_W = 0;
mix_theta = 0;

%% Object containing EM optimization optimization
EM = EMstats;
EM = EM.setMaxIterations(30);
EM = EM.prealloc(fineData, nFine, nCoarse, nBasis);           %preallocation of data arrays





