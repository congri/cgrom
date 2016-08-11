%Main script for coarse-graining/reduced order modeling
clear all;
addpath('./MCMCsampler')
addpath('./util')

%random number seed
rng(0)

tic;
%load params
params;

Fmesh = genMesh(boundary, nFine);
Cmesh = genMesh(boundary, nCoarse);

%Generate finescale dataset
[x, Tf, PhiArray] = genFineData(Fmesh, phi, heatSource, boundary, fineCond, nFine, nCoarse);

% If no parallel pool exists
N_Threads = 5;
if isempty(gcp('nocreate'))
    % Create with N_Threads workers
    parpool('local',N_Threads);
end

%store handle to every q_i in a cell array lq
lq = cell(fineCond.nSamples, 1);

sw = zeros(fineCond.nSamples, maxIterations);
for k = 1:maxIterations
    %Generate samples from every q_i
    parfor i = 1:fineCond.nSamples
        lq{i} = @(Xi) log_q_i(Xi, Tf(:,i), theta_cf, theta_c, PhiArray(:,:,i), Fmesh, Cmesh, heatSource, boundary, W);
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
        
        %Compute sufficient statistics
        if any(any(out(i).samples < 0))
            error('negative components')
        end
        XMean(:, i) = mean(out(i).samples, 2);
        XNormSqMean(i) = mean(sum(out(i).samples.^2));
        temp(:,i) = cell2mat(out(i).data);
        muCfSq(:,:,i) = reshape(temp(:,i).^2, nFine + 1, MCMC(i).nSamples);
        muCfSqMean(:,i) = mean(muCfSq(:,:,i), 2);
    end
    lp
    theta_cf.S = diag(mean(muCfSqMean, 2));
    %ensure invertability; noise vanishes at essential nodes
    stabilityFactor = 1e-9;
    if(~theta_cf.S(1))
        theta_cf.S(1) = stabilityFactor;
    end
    if(~theta_cf.S(end))
        theta_cf.S(end) = stabilityFactor;
    end
    
    sumPhiSq = zeros(size(phi, 1), size(phi, 1));
    sumPhiTXmean = zeros(size(phi, 1), 1);
    for i = 1:fineCond.nSamples
        sumPhiSq = sumPhiSq + PhiArray(:,:,i)'*PhiArray(:,:,i);
        sumPhiTXmean = sumPhiTXmean + PhiArray(:,:,i)'*XMean(:,i);
    end
    theta_c.theta = sumPhiSq\sumPhiTXmean;
    
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
    theta_c.sigma = sqrt(sigmaSq);
    thetaArray(:,k) = theta_c.theta
    sigmaArray(k) = theta_c.sigma
    S = diag(theta_cf.S)
    
end


runtime = toc




