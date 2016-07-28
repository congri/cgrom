function [F] = getForce(domain, heatSource, boundary)
%Compute the global force vector

%There are only equations for natural nodes
F = domain.nodes(domain.nodes == 0);

%Heat source part
if strcmp(heatSource.type, 'const')

    %heat source
    F(1) = F(1) + .5*domain.l*heatSource.value;
    F(end) = F(end) + .5*domain.l*heatSource.value;
    
    %essential boundary
    F(1) = F(1) - boundary.T0(1)*domain.nodes(1)*domain.conductivity(1)/domain.l;
    F(2) = F(2) + boundary.T0(1)*domain.nodes(1)*domain.conductivity(1)/domain.l;
    F(end - 1) = F(end - 1) + boundary.T0(2)*domain.nodes(end)*domain.conductivity(end)/domain.l;
    F(end) = F(end) - boundary.T0(2)*domain.nodes(end)*domain.conductivity(end)/domain.l;
    
    nF = numel(F);
    if nF > 2
        for i = 2:(nF - 1)
            F(i) = F(i) + domain.l*heatSource.value;
        end
    end

end



end

