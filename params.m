%parameters for 1D FEM
domain.N_el = 3;                              %number of elements
domain.l = 1/domain.N_el;                             %element length
domain.conductivity = linspace(1,10,domain.N_el);     %conductivities

heatSource.type = 'const';               %type of heat source field; const is const in space
heatSource.value = 1;                                  %heat source field

%boundary conditions, 1 element is left, 2 is right boundary
boundary.type = {'essential', 'natural'};
boundary.T0 = [5; 0];
boundary.q0 = [0; 1];

%Essential (1) or natural (0) nodes?
domain.nodes = zeros(domain.N_el + 1, 1);
if strcmp(boundary.type(1), 'essential')
    domain.nodes(1) = 1;
end
if strcmp(boundary.type(2), 'essential')
    domain.nodes(end) = 1;
end


