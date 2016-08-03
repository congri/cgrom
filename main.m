%Main script for coarse-graining/reduced order modeling

%get params
domain.N_el = 4;
params;

%compute global force vector
F = getForce(domain, heatSource, boundary);

K = getStiff(domain);

%nodal temperatures
Tf = K\F;
if domain.nodes(1)
    Tf = [boundary.T0(1); Tf];
end
if domain.nodes(end)
    Tf = [Tf; boundary.T0(2)];
end
Tf


%get params
domain.N_el = 2;
params;

%compute global force vector
F = getForce(domain, heatSource, boundary);

K = getStiff(domain);

%nodal temperatures
Tc = K\F;
if domain.nodes(1)
    Tc = [boundary.T0(1); Tc];
end
if domain.nodes(end)
    Tc = [Tc; boundary.T0(2)];
end
Tc