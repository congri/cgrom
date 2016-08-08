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
[x, Tf, PhiArray] = genFineData(Fmesh, phi, heatSource, boundary, fineCond, nFine, nCoarse);

% If no parallel pool exists
N_Threads = 2;
if isempty(gcp('nocreate'))
    % Create with 2 workers
    parpool('local',N_Threads);
end

%store handle to every q_i in a cell array lq
lq = cell(fineCond.nSamples, 1);

for k = 1:100
    %Generate samples from every q_i
    parfor i = 1:fineCond.nSamples
        lq{i} = @(Xi) log_q_i(Xi, Tf(:,i), theta_cf, theta_c, PhiArray(:,:,i), Fmesh, Cmesh, heatSource, boundary, W);
        %sample from every q_i
        out(i) = MCMCsampler(lq{i}, MCMC(i).Xi_start, MCMC(i));
        out(i).log_pEnd
        MCMC(i).Xi_start = out(i).samples(:, end);
        
        %Refine step width
        if(out(i).acceptance)
            MCMC(i).randomWalk.proposalCov = (1/.6)*out(i).acceptance*MCMC(i).randomWalk.proposalCov;
        else
            MCMC(i).randomWalk.proposalCov = .2*MCMC(i).randomWalk.proposalCov;
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
    
    sigmaSq = 0;
    for i = 1:fineCond.nSamples
        sigmaSq = sigmaSq + XNormSqMean(i) - 2*theta_c.theta'*PhiArray(:,:,i)'*XMean(:,i)...
            + theta_c.theta'*PhiArray(:,:,i)'*PhiArray(:,:,i)*theta_c.theta;
    end
    sigmaSq = sigmaSq/(nCoarse*fineCond.nSamples);
    theta_c.sigma = sqrt(sigmaSq);
    thetaArray(:,k) = theta_c.theta
    sigmaArray(k) = theta_c.sigma
    S = theta_cf.S
    
end


runtime = toc




