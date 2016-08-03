%parameters for 1D FEM
heatSource.type = 'const';                              %type of heat source field; const is const in space
heatSource.value = 0;                                   %heat source field

%boundary conditions, 1 element is left, 2 is right boundary
boundary.type = {'essential', 'natural'};
boundary.T0 = [0; 0];
boundary.q0 = [0; 1];

%Finescale conductivity params
fineCond.mu = 1;    %mean of log of lambda
fineCond.sigma = 1; %sigma of log of lambda
fineCond.Nsamples = 10;

%Fine and coarse number of elements
nFine = 9;
nCoarse = 3;
assert(~mod(nFine, nCoarse), 'Error: Coarse mesh is not a divisor of fine mesh!')

%Define basis functions for p_c here
phi_1 = @(x) 1/sum(1./x);
phi_2 = @(x) mean(x);
phi = {phi_1; phi_2};



