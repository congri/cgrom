function [T] = FEMmain(domain, heatSource, boundary)
%Core of FEM simulation

%compute global force vector
[F, d_F] = getForce(domain, heatSource, boundary);

[K, d_K] = getStiff(domain);

%nodal temperatures
T = K\F;
if domain.nodes(1)
    T = [boundary.T0(1); T];
end
if domain.nodes(end)
    T = [T; boundary.T0(2)];
end
    
    
end

