%% Main script for coarse-graining/reduced order modeling
%% preamble
clear all;
addpath('./MCMCsampler')
addpath('./util')
addpath('./FEM')
addpath('./model')
addpath('./params')

%random number seed
% rng(0)
rng('shuffle')
tic;

%% load params
params;

%% Generate fine and coarse meshes
Fmesh = genMesh(boundary, nFine);
Cmesh = genMesh(boundary, nCoarse);

%% Generate finescale dataset
if fineCond.genData
    [x, Tf, PhiArray] = genFineData(Fmesh, phi, heatSource, boundary, fineCond, nFine, nCoarse);
else
    load('./data/fineData/fineData')
end
sumPhiSq = zeros(size(phi, 1), size(phi, 1));
for i = 1:fineCond.nSamples
    sumPhiSq = sumPhiSq + PhiArray(:,:,i)'*PhiArray(:,:,i);
    %take MCMC initializations at mode of p_c
    MCMC(i).Xi_start = PhiArray(:,:,i)*theta_c.theta;
end
sumPhiSqInv = inv(sumPhiSq);

%% Create parallel pool
N_Threads = min(15, fineCond.nSamples);
if isempty(gcp('nocreate'))
    % Create with N_Threads workers
    parpool('local',N_Threads);
end

%% EM optimization - main body
%store handle to every q_i in a cell array lq
lq = cell(fineCond.nSamples, 1);

sw = zeros(fineCond.nSamples, maxIterations);
Warray = 0*repmat(theta_cf.W,1,1,maxIterations);
for k = 1:maxIterations
    
    Warray(:,:,k ) = theta_cf.W;
    Wplt(:,k) = theta_cf.W(:);
    
    %Generate samples from every q_i
    parfor i = 1:fineCond.nSamples
        lq{i} = @(Xi) log_q_i(Xi, Tf(:,i), theta_cf, theta_c, PhiArray(:,:,i), Fmesh, Cmesh, heatSource, boundary, theta_cf.W);
        %sample from every q_i
        out(i) = MCMCsampler(lq{i}, MCMC(i).Xi_start, MCMC(i));
        %avoid very low acceptances
        while out(i).acceptance < .1
            out(i) = MCMCsampler(lq{i}, MCMC(i).Xi_start, MCMC(i));
            if strcmp(MCMC(i).method, 'MALA')
                MCMC(i).MALA.stepWidth = (1/.8)*(out(i).acceptance + (1 - out(i).acceptance)*.1)*MCMC(i).MALA.stepWidth;
            elseif strcmp(MCMC(i).method, 'randomWalk')
                MCMC(i).randomWalk.proposalCov = .8*MCMC(i).randomWalk.proposalCov;
            else
            end
            warning('Acceptance ratio below .1')
        end
        sw(i, k) = MCMC(i).MALA.stepWidth;
        %         MCMC(i).Xi_start = out(i).samples(:, end);
        
        lp(i) = mean(out(i).log_p);
        
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
        %temp contains coarse temperatures for all coarse nodes and all fine data samples
        temp(:,i) = cell2mat(out(i).data);
        %Tc_samples(:,:,i) contains coarse nodal temperature samples (1 sample == 1 column) for data
        %sample i
        Tc_samples(:,:,i) = reshape(temp(:,i), nCoarse + 1, MCMC(i).nSamples);
        %only valid for diagonal S here!
        temp2(:, i) = mean((repmat(Tf(:, i) - theta_cf.mu, 1, MCMC(i).nSamples) - theta_cf.W*Tc_samples(:,:,i)).^2, 2);
        
        %for W
        %sample mean under q_i
        Tc_mean(:,i) = mean(Tc_samples(:,:,i), 2);
        %First factor for matrix W
        Wa(:,:,i) = (Tf(:, i) - theta_cf.mu)*Tc_mean(:, i)';
    end
    
    for i = 1:fineCond.nSamples
        for l = 1:MCMC(i).nSamples
            Tc_dyadic(:,:,l) = Tc_samples(:,l,i)*Tc_samples(:,l,i)';
        end
        %mean of dyadic product under q_i
        Tc_dyadic_mean(:,:,i) = mean(Tc_dyadic, 3);
    end
    %Mean along data samples
    Tc_dyadic_mean_mean = mean(Tc_dyadic_mean, 3);
    Wa_mean = mean(Wa, 3);
    if(~Winterp)
        %     theta_cf.W = (1 - mix_W)*(Wa_mean/Tc_dyadic_mean_mean) + mix_W*theta_cf.W;
        theta_cf.W = compW(Tc_dyadic_mean_mean, Wa_mean,...
            inv(theta_cf.S), theta_cf.W, paramIndices, constIndices);
    end
    
%     Wout = theta_cf.W
    %     lp
    %decelerate convergence of S
    theta_cf.S = (1 - mix_S)*diag(mean(temp2, 2)) + mix_S*theta_cf.S + (1e-6)*eye(size(theta_cf.S, 1));
    
    
    sumPhiTXmean = zeros(size(phi, 1), 1);
    for i = 1:fineCond.nSamples
        sumPhiTXmean = sumPhiTXmean + PhiArray(:,:,i)'*XMean(:,i);
    end
    theta_c.theta = (1 - mix_theta)*(sumPhiSqInv*sumPhiTXmean) + mix_theta*theta_c.theta;
    
    %Start next chain at mean of p_c
    for i = 1:fineCond.nSamples
        MCMC(i).Xi_start = PhiArray(:,:,i)*theta_c.theta;
    end
    
    sigmaSq = 0;
    for i = 1:fineCond.nSamples
        sigmaSq = sigmaSq + XNormSqMean(i) - 2*theta_c.theta'*PhiArray(:,:,i)'*XMean(:,i)...
            + theta_c.theta'*PhiArray(:,:,i)'*PhiArray(:,:,i)*theta_c.theta;
    end
    sigmaSq = sigmaSq/(nCoarse*fineCond.nSamples);
    sigmaOffset = 1e-8;
    theta_c.sigma = (1 - mix_sigma)*sqrt(sigmaSq) + mix_sigma*theta_c.sigma + sigmaOffset;
    thetaArray(:,k) = theta_c.theta
    sigmaArray(k) = theta_c.sigma
    S = diag(theta_cf.S)'
    
end


runtime = toc




