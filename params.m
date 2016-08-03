%parameters for 1D FEM
%domain.N_el = 100;                                        %number of elements
domain.l = 1/domain.N_el;                                  %element length
heatSource.type = 'const';                              %type of heat source field; const is const in space
heatSource.value = 0;                                   %heat source field

%boundary conditions, 1 element is left, 2 is right boundary
boundary.type = {'essential', 'natural'};
boundary.T0 = [0; 0];
boundary.q0 = [0; 1];

%Essential (1) or natural (0) nodes?
domain.nodes = false(domain.N_el + 1, 1);
if strcmp(boundary.type(1), 'essential')
    domain.nodes(1) = true;
end
if strcmp(boundary.type(2), 'essential')
    domain.nodes(end) = true;
end

%Finescale conductivity params
fineCond.mu = 1;    %mean of log of lambda
fineCond.sigma = 1; %sigma of log of lambda
fineCond.Nsamples = 10;


