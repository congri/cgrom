function [curr_x, grad, Hess, nIter] = newtonRaphsonMaximization(objective, startValue, Xtol, provide_objective)
%Maximization based on the well-known Newton-Raphson method
%%Input:
%   objective:              func. handle with [grad., Hess., obj.] as output
%   startValue:             Optimization start value
%   Xtol:                   Tolerance in x for convergence

%%Optimization
converged = false;
curr_x = startValue;        %clear notation
if provide_objective
    [grad, Hess, obj] = objective(curr_x);
else
    [grad, Hess] = objective(curr_x);
end
step = Hess\grad;
nIter = 1;
%Iterate until convergence
while(~converged)
    old_x = curr_x;
    if(~provide_objective)
        curr_x = old_x - step;
    else
        temp_x = old_x - step;
        [grad_temp, Hess_temp, obj_temp] = objective(temp_x)
        nIter = nIter + 1;
        step_temp = Hess_temp\grad_temp;
        while(obj_temp < obj)
            step_temp = .9*step_temp;
            temp_x = old_x - step_temp;
            [grad_temp, Hess_temp, obj_temp] = objective(temp_x);
            nIter = nIter + 1;
        end
        %step accepted
        curr_x = temp_x;
        obj = obj_temp;
        grad = grad_temp;
        Hess = Hess_temp;
        step = Hess\grad;
    end
    
    %check for convergence
    if(norm(curr_x - old_x)/norm(curr_x) < Xtol)
        converged = true;
    else    %if not converged, compute next step
        if(~provide_objective)
            [grad, Hess] = objective(curr_x);
            nIter = nIter + 1;
            step = Hess\grad;
        end
    end
end

end

