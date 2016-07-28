%parameters for 1D FEM
domain.N_el = 3;                              %number of elements
domain.l = 1/domain.N_el;                             %element length
domain.conductivity = linspace(1,10,domain.N_el);     %conductivities

heatSource.type = 'const';               %type of heat source field; const is const in space
heatSource.value = 0;                                  %heat source field

%boundary conditions, 1 element is left, 2 is right boundary
boundary.T0 = [5; NaN];
boundary.q0 = [NaN; 1];

%Essential (1) or natural (0) nodes?
domain.nodes = zeros(domain.N_el + 1, 1);
essBoundaries = isnan(boundary.q0);
if essBoundaries(1)
    domain.nodes(1) = 1;
end
if essBoundaries(2)
    domain.nodes(end) = 1;
end


