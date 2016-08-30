%% Main script for coarse-graining/reduced order modeling
%% preamble
clear all;
addpath('./MCMCsampler')
addpath('./util')
addpath('./FEM')
addpath('./model')
addpath('./params')
addpath('./computation')

rng('shuffle')
tic;

%% load params
params;

%% Generate fine and coarse meshes
Fmesh = genMesh(boundary, nFine);
Cmesh = genMesh(boundary, nCoarse);

%% Generate finescale dataset
if fineData.genData
    [x, Tf, PhiArray, sumPhiSqInv] = genFineData(Fmesh, phi, heatSource, boundary, fineData, nFine, nCoarse);
else
    load('./data/fineData/fineData')
end
for i = 1:fineData.nSamples
    %take MCMC initializations at mode of p_c
    MCMC(i).Xi_start = PhiArray(:,:,i)*theta_c.theta;
end

%% Open parallel pool
parPoolInit(fineData.nSamples);

%% EM optimization - main body
%store handle to every q_i in a cell array lq
log_qi = cell(fineData.nSamples, 1);

%collect data in data arrays
k = 1;
collectData;
for k = 2:(EM.maxIterations + 1)    
    %Generate samples from every q_i
    parfor i = 1:fineData.nSamples
        log_qi{i} = @(Xi) log_q_i(Xi, Tf(:,i), theta_cf, theta_c, PhiArray(:,:,i), Fmesh, Cmesh, heatSource, boundary, theta_cf.W);
        %sample from every q_i
        out(i) = MCMCsampler(log_qi{i}, MCMC(i).Xi_start, MCMC(i));
        %avoid very low acceptances
        while out(i).acceptance < .1
            out(i) = MCMCsampler(log_qi{i}, MCMC(i).Xi_start, MCMC(i));
            if strcmp(MCMC(i).method, 'MALA')
                MCMC(i).MALA.stepWidth = (1/.8)*(out(i).acceptance + (1 - out(i).acceptance)*.1)*MCMC(i).MALA.stepWidth;
            elseif strcmp(MCMC(i).method, 'randomWalk')
                MCMC(i).randomWalk.proposalCov = .8*MCMC(i).randomWalk.proposalCov;
            else
                error('Unknown MCMC method')
            end
            warning('Acceptance ratio below .1')
        end
        
%         log_qi_mean(i) = mean(out(i).log_p);
        
        %Refine step width
        if strcmp(MCMC(i).method, 'MALA')
            MCMC(i).MALA.stepWidth = (1/.8)*out(i).acceptance*MCMC(i).MALA.stepWidth;
        elseif strcmp(MCMC(i).method, 'randomWalk')
            MCMC(i).randomWalk.proposalCov = (1/.5)*out(i).acceptance*MCMC(i).randomWalk.proposalCov;
        else
        end
        
        XMean(:, i) = mean(out(i).samples, 2);
        XNormSqMean(i) = mean(sum(out(i).samples.^2));
        
        %for S
        %Tc_samples(:,:,i) contains coarse nodal temperature samples (1 sample == 1 column) for full order data
        %sample i
        Tc_samples(:, :, i) = reshape(cell2mat(out(i).data), nCoarse + 1, MCMC(i).nSamples);
        %only valid for diagonal S here!
        p_cf_exponent(:, i) = mean((repmat(Tf(:, i) - theta_cf.mu, 1, MCMC(i).nSamples) - theta_cf.W*Tc_samples(:, :, i)).^2, 2);
        
        %First factor for matrix W
        Wa(:, :, i) = (Tf(:, i) - theta_cf.mu)*mean(Tc_samples(:, :, i), 2)';
    end
    
    Tc_dyadic_mean = TcDyadicMean(Tc_samples, fineData.nSamples, MCMC);
    Wa_mean = mean(Wa, 3);
    if(~Winterp)
        theta_cf.W = compW(Tc_dyadic_mean,  Wa_mean,...
            inv(theta_cf.S), theta_cf.W, paramIndices, constIndices);
    end
    
%     Wout = theta_cf.W
    %decelerate convergence of S
    theta_cf.S = (1 - mix_S)*diag(mean(p_cf_exponent, 2)) + mix_S*theta_cf.S + (1e-6)*eye(size(theta_cf.S, 1));
    
    
    sumPhiTXmean = zeros(size(phi, 1), 1);
    for i = 1:fineData.nSamples
        sumPhiTXmean = sumPhiTXmean + PhiArray(:,:,i)'*XMean(:,i);
    end
    theta_c.theta = (1 - mix_theta)*(sumPhiSqInv*sumPhiTXmean) + mix_theta*theta_c.theta;
    curr_theta = theta_c.theta
    
    %Start next chain at mean of p_c
    for i = 1:fineData.nSamples
        MCMC(i).Xi_start = PhiArray(:,:,i)*theta_c.theta;
%         MCMC(i).Xi_start = out(i).samples(:, end);
    end
    
    sigmaSq = 0;
    for i = 1:fineData.nSamples
        sigmaSq = sigmaSq + XNormSqMean(i) - 2*theta_c.theta'*PhiArray(:,:,i)'*XMean(:,i)...
            + theta_c.theta'*PhiArray(:,:,i)'*PhiArray(:,:,i)*theta_c.theta;
    end
    sigmaSq = sigmaSq/(nCoarse*fineData.nSamples);
    sigmaOffset = 1e-8;
    theta_c.sigma = (1 - mix_sigma)*sqrt(sigmaSq) + mix_sigma*theta_c.sigma + sigmaOffset;
    
    S = diag(theta_cf.S)'
    %% collect data in data arrays
    collectData;
end
clear i j k m Wa Wa_mean Tc_dyadic_mean log_qi out p_cf_exponent curr_theta XMean XNormSqMean;


runtime = toc




