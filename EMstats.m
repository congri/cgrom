classdef EMstats
    %Class describing the Expectation-Maximization procedure and storing the data sequentially
    
    properties
        %% Optimization data
        MCMCStepWidth
        
        %% Optimal params, theta_cf
        W
        S
        mu
        %% Optimal params, theta_c
        theta
        sigma
        
    end
    
    properties (SetAccess = private)
        %read-only properties
        
        %Optimization params
        maxIterations
        
    end
    
    methods
        
        function obj = setMaxIterations(obj, maxIterations)
            %% some tests
            assert((maxIterations > 0), 'Number of EM iterations has to be > 0')
            assert((mod(maxIterations, 1) == 0), 'Number of EM iterations must be integer')
            assert(all(size(maxIterations)) == 1, 'maximum number of iterations must be a scalar')
            
            %% set value
            obj.maxIterations = int32(maxIterations);
        end
        
        function obj = prealloc(obj, fineData, nFine, nCoarse, nBasis)
            %preallocation of data arrays
            
            %% tests
            assert(~mod(nFine, nCoarse), 'Error: Coarse mesh is not a divisor of fine mesh!')
            
            %% for data record, single precision is sufficient
            obj.MCMCStepWidth = zeros(fineData.nSamples, obj.maxIterations + 1, 'single');
            obj.W = zeros(nFine + 1, nCoarse + 1, obj.maxIterations, 'single');
            obj.theta = zeros(nBasis, obj.maxIterations + 1, 'single');
            obj.sigma = zeros(obj.maxIterations + 1, 1, 'single');
            obj.S = zeros(nFine + 1, nFine + 1, obj.maxIterations + 1, 'single');
            obj.mu = zeros(nFine + 1, obj.maxIterations, 'single');
        end
        
        function plotData(obj)
            %plot the final data
            f = figure;
            set(f, 'Position', [680 678 900 700]) %[position x, psition y, width, height]
            subplot(2, 2, 1)
            iterations = 1:(obj.maxIterations + 1);
            p_theta = plot(iterations, obj.theta);
            set(gca, 'fontsize', 15)
            set(p_theta, 'linewidth', 2)
            xlabel('iteration i')
            ylabel('\theta_c')
            xlim([1 obj.maxIterations + 1])
            ylim([-3 3])
            axis square
            grid on
            
            subplot(2, 2, 2)
            p_sigma = plot(iterations, obj.sigma);
            set(gca, 'fontsize', 15)
            set(p_sigma, 'linewidth', 2)
            axis square
            grid on
            xlabel('iteration i')
            ylabel('\sigma')
            xlim([1 obj.maxIterations + 1])
            
            
            subplot(2, 2, 3)
            Splt = zeros(size(obj.S, 1), obj.maxIterations + 1);
            for i = 1:(obj.maxIterations + 1)
                Splt(:, i) = diag(obj.S(:, :, i));
            end
            p_S = plot(iterations, Splt, 'k');
            set(gca, 'fontsize', 15)
            set(p_S, 'linewidth', 1)
            axis square
            grid on
            xlabel('iteration i')
            ylabel('S')
            xlim([1 obj.maxIterations + 1])
            
            subplot(2, 2, 4)
            
            p_sw = semilogy(iterations, obj.MCMCStepWidth, 'k');
            set(gca, 'fontsize', 15)
            set(p_sw, 'linewidth', 1)
            axis square
            grid on
            xlabel('iteration i')
            ylabel('step width')
            xlim([1 obj.maxIterations + 1])
            
        end
        
    end
    
end
