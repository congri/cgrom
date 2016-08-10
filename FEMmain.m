function [T, d_r, K] = FEMmain(domain, heatSource, boundary)
%Core of FEM simulation
%d_r is the derivative of the whole equation system

%compute global force vector
[F, d_F] = getForce(domain, heatSource, boundary);
[K, d_K] = getStiff(domain);

%nodal temperatures
T = K\F;

%Compute derivative of FEM system
d_r = zeros(size(F, 1), domain.N_el);
for i = 1:domain.N_el  
   d_r(:,i) = d_K(:,:,i)*T - d_F(:,i);
end

if domain.nodes(1)
    T = [boundary.T0(1); T];
end
if domain.nodes(end)
    T = [T; boundary.T0(2)];
end
    
    
end

