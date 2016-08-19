function [F, d_F] = getForce(domain, heatSource, boundary)
%Compute the global force vector

%There are only equations for natural nodes
F = double(domain.nodes(domain.nodes == 0));
d_F = repmat(F, 1, domain.N_el);

%Heat source part
if strcmp(heatSource.type, 'const')

    %heat source
    F(1) = F(1) + .5*domain.l*heatSource.value;
    F(end) = F(end) + .5*domain.l*heatSource.value;
    
    nF = numel(F);
    if nF > 2
        for i = 2:(nF - 1)
            F(i) = F(i) + domain.l*heatSource.value;
        end
    end
    
    %boundaries
    if domain.nodes(1)
        %left end is essential
        %keep in mind that 1st equation corresponds to 2nd node
        F(1) = F(1) + domain.conductivity(1)*boundary.T0(1)/domain.l;
        d_F(1, 1) = d_F(1, 1) + boundary.T0(1)/domain.l;
    else
        %left end is natural
        F(1) = F(1) + boundary.q0(1);
    end
    if domain.nodes(end)
        %right end is essential
        %keep in mind that last equation corresponds to next to last node
        F(end) = F(end) + domain.conductivity(end)*boundary.T0(2)/domain.l;
        d_F(end, end) = d_F(end, end) + boundary.T0(2)/domain.l;
    else
        %right end is natural
        F(end) = F(end) - boundary.q0(2);
    end

    
else
    
    error('Unknown heat source type')
    
end


end

