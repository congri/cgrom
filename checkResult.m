function [] = checkResult(theta_cIn, theta_cfIn, phiIn)
%We sample the resulting distribution p(y|x, theta_c, theta_cf) here and compare to the true
%solution

%read params
params;
clear MCMC;
theta_c = theta_cIn;
theta_cf = theta_cfIn;
phi = phiIn;
%MCMC options
MCMC.method = 'MALA';                             %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 100;
MCMC.nThermalization = 3000;                              %thermalization steps
MCMC.nSamples = 20;                                    %number of samples
MCMC.nGap = 5000;
MCMC.Xi_start = ones(nCoarse, 1);
%only for random walk
MCMC.randomWalk.proposalCov = stepWidth*eye(nCoarse);   %random walk proposal covariance

fineData.nSamples = 1;
%type optimal parameters here
MCMC.MALA.stepWidth = theta_c.sigma;

%generate reference data
boundary.type;
Fmesh = genMesh(boundary, nFine);
%rng(1)
[x, Tf] = genFineData(Fmesh, heatSource, boundary, fineData);
Phi = designMatrix(phi, exp(x), nFine, nCoarse);

%samples from p_c
l_p_c = @(X) log_p_c(X, Phi, theta_c.theta, theta_c.sigma);
out_p_c = MCMCsampler(l_p_c, MCMC.Xi_start, MCMC);

%solve CG model
MCMC_cf = MCMC;
MCMC_cf.Xi_start = Tf(:,1);
MCMC_cf.method = 'randomWalk';
MCMC_cf.nThermalization = 3000;
MCMC_cf.nSamples = 1;
stepWidth = 1e-1;
MCMC_cf.randomWalk.proposalCov = theta_cf.S;   %random walk proposal covariance

Cmesh = genMesh(boundary, nCoarse);
for i = 1:MCMC.nSamples
    i
    Cmesh.conductivity = out_p_c.samples(:, i);
    
    %draw from p_cf
    Tf
    l_p_cf = @(Tfv) log_p_cf(Tfv, Cmesh, heatSource, boundary, theta_cf.W, theta_cf.S);
    out_p_cf = MCMCsampler(l_p_cf, MCMC_cf.Xi_start, MCMC_cf);
    out_p_cf
    out_p_cf.samples
    if i == 1
        y_samples = out_p_cf.samples;
    else
        y_samples = [y_samples, out_p_cf.samples];
    end
end

figure
pt = plot(Tf);
hold
ps = plot(y_samples(:,1), 'xk');
legend('FOM solution', 'model samples')
plot(y_samples(:,2:end), 'xk')
axis square
grid on
xlabel('node i')
ylabel('temperature y_i')







