% Define your objective function
fitnessFunction = @(x) sum(x.^2); % Example objective function, replace with your own

% Define the number of variables
nvars = 5; % Example: 5 variables

% Define the options for the genetic algorithm
options = optimoptions('ga', 'OutputFcn', @outputFunction, 'Display', 'iter');

% Call the ga function with your objective function, number of variables, lower and upper bounds, and options
[x, fval] = ga(fitnessFunction, nvars, [], [], [], [], zeros(nvars, 1), ones(nvars, 1), [], options);

% Define the custom output function to track progress
function [state,options,optchanged] = outputFunction(options,state,flag)
    persistent bestFitnessHistory;
    switch flag
        case 'init'
            bestFitnessHistory = [];
        case 'iter'
            % Store the best fitness value at each generation
            bestFitnessHistory = [bestFitnessHistory, min(state.Score)];
            % Display the best fitness value at each generation
            disp(['Generation ', num2str(state.Generation), ': Best Fitness = ', num2str(min(state.Score))]);
        case 'done'
            % Display final message
            disp('GA optimization complete.');
    end
    optchanged = false; % Continue optimization
end
