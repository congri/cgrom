%test script for probability distribution q_i

%random number seed
rng(1)

%boundary conditions, 1 element is left, 2 is right boundary
boundary.type = {'essential', 'natural'};
boundary.T0 = [3; 1];
boundary.q0 = [1; 1];

heatSource.type = 'const';                              %type of heat source field; const is const in space
heatSource.value = 0; 

nFine = 4;
nCoarse = 2;

Fmesh = genMesh(boundary, nFine);
Cmesh = genMesh(boundary, nCoarse);

%Define basis functions for p_c here
phi_1 = @(x) size(x, 1)/sum(1./x);
phi_2 = @(x) mean(x);
phi = {phi_1};

%Finescale conductivity params
fineCond.mu = 1;    %mean of log of lambda
fineCond.sigma = 1; %sigma of log of lambda
fineCond.nSamples = 1;

%fine data
[x, Tf, Phi] = genFineData(Fmesh, phi, heatSource, boundary, fineCond, nFine, nCoarse);

%start values
theta_cf.S = 1*eye(nFine + 1);
theta_c.theta = [1];
theta_c.sigma = 10;

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


x = linspace(.1, 40, 50);
[X1, X2] = meshgrid(x);
Xi = [X1(:)'; X2(:)'];

lq = 0*X1;
for i = 1:size(Xi, 2)
    [log_q, d_log_q, mu_pcf] = log_q_i(Xi(:,i), Tf, theta_cf, theta_c, Phi, Fmesh, Cmesh, heatSource, boundary, W);
    lq(i) = log_q;
end

figure
subplot(1,2,1)
contourf(X1, X2, lq)
axis square
subplot(1,2,2)
contourf(X1, X2, exp(lq))
axis square
