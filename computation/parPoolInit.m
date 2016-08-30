function [] = parPoolInit(nDataSamples)
%% Initializes parallel pool

nLocalCores = 16;   %number of cores on local machine
N_Threads = min(nLocalCores, nDataSamples);
if isempty(gcp('nocreate'))
    % Create with N_Threads workers
    parpool('local', N_Threads);
end

end

