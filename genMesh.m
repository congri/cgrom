function [domain] = genMesh(boundary, N_el)
%Generates a FE mesh with the given number of elements

domain.N_el = N_el;
domain.l = 1/N_el;                                  %element length

%Essential (1) or natural (0) nodes?
domain.nodes = false(domain.N_el + 1, 1);
if strcmp(boundary.type(1), 'essential')
    domain.nodes(1) = true;
end
if strcmp(boundary.type(2), 'essential')
    domain.nodes(end) = true;
end


    
end

