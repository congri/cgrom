function [phi] = secOrderFeature(x)
%second order feature function

phi = 0;
for i = 1:(length(x) - 1)
    phi = phi + x(i)*x(i + 1);   
end


end

