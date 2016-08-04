
N_Threads = 2;

N = 500;
x = randn(N,1);
y = zeros(N,1);
dt = 0.01;

%% SERIAL

disp('   Run in Serial');

tic
for i=1:N
    % call function
    y(i) = sin(x(i));
    % pretend like we're doing something remotely computational expensive
    % (otherwise parallelization overhead would dominate);
    pause(dt);
end
t_s = toc;


disp('   Run in Parallel');

% If no parallel pool exists
if isempty(gcp('nocreate'))
    % Create with 2 workers
    parpool('local',N_Threads);
end

%% PARALLEL

tic
parfor i=1:N
    % call function
    y(i) = sin(x(i));
    % pretend like we're doing something remotely computational expensive    
    pause(dt);
end
t_p = toc;

disp(['Runtime Serial: ' num2str(t_s)]);
disp(['Runtime Parallel: ' num2str(t_p)]);
disp(['Speedup: ' num2str(t_s/t_p)]);

gcp



