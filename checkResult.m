%We sample the resulting distribution p(y|x, theta_c, theta_cf) here and compare to the true
%solution

%Main script for coarse-graining/reduced order modeling
clear all;
addpath('./MCMCsampler')
addpath('./util')
addpath('./FEM')
addpath('./model')
addpath('./params')


%read params
params;
clear MCMC;
%MCMC options
MCMC.method = 'MALA';                             %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 100;
MCMC.nThermalization = 3000;                              %thermalization steps
MCMC.nSamples = 20;                                    %number of samples
MCMC.nGap = 5000;
MCMC.Xi_start = ones(nCoarse, 1);
%only for random walk
MCMC.randomWalk.proposalCov = stepWidth*eye(nCoarse);   %random walk proposal covariance

fineCond.nSamples = 1;
%type optimal parameters here
theta_c.theta = [.1658
                .00066372];
theta_c.sigma = .021801;
MCMC.MALA.stepWidth = theta_c.sigma;


theta_cf.W = [1   2.4404e-16   1.2156e-15  -1.0733e-15
      0.80365      0.10814    -0.022919    0.0019982
      0.83083      0.16551   0.00021975  -0.00039017
      0.80476      0.10361     0.047167     0.011653
      0.82775  -0.00052949       0.1674   -0.0012295
       1.0218     0.033287     0.036077      0.12347
      0.83394   -0.0001552    -0.004427       0.1698];
  
theta_cf.S = [1.3287e-28            0            0            0            0            0            0
            0       1.0153            0            0            0            0            0
            0            0    2.284e-05            0            0            0            0
            0            0            0       2.6667            0            0            0
            0            0            0            0   0.00018717            0            0
            0            0            0            0            0       2.0741            0
            0            0            0            0            0            0   8.9364e-06];

%generate reference data
Fmesh = genMesh(boundary, nFine);
%rng(1)
[x, Tf, Phi] = genFineData(Fmesh, phi, heatSource, boundary, fineCond, nFine, nCoarse);

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
    l_p_cf = @(Tfv) log_p_cf(Tfv, Cmesh, heatSource, boundary, theta_cf.W, theta_cf.S);
    out_p_cf = MCMCsampler(l_p_cf, MCMC_cf.Xi_start, MCMC_cf);
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







