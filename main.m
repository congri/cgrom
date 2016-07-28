%Main script for coarse-graining/reduced order modeling

%get params
params;

%compute global force vector
F = getForce(domain, heatSource, boundary)