% Define data for each item
weights = [5, 1, 9, 9, 6, 5, 2, 5, 6, 7];
volumes = [12, 12, 17, 14, 15, 12, 19, 12, 20, 13];
values = [62.30, 32.53, 59.38, 65.12, 48.67, 14.36, 69.52, 84.71, 6.58, 46.99];
nItems = length(weights);

% Constraints
maxWeightLimit = 30;
maxVolumeLimit = 60;
maxItemsSelected = 5;

% Variable bounds (binary representation of item selection)
lb = zeros(nItems, 1); % Lower bound
ub = ones(nItems, 1);  % Upper bound

% Integer constraint for GA (all variables are integers)
IntCon = 1:nItems;

% Objective function for the knapsack problem
objfun = @(x) -sum(values .* x); % Maximize the total value

% Constraint function for the knapsack problem
confun = @(x) knapsackConstraints(x, weights, volumes, maxWeightLimit, maxVolumeLimit, maxItemsSelected);

% GA options
options = optimoptions('ga', ...
    'PopulationSize', 200, ...
    'MaxGenerations', 150, ...
    'CrossoverFraction', 0.8, ...
    'MutationFcn', {@mutationuniform, 0.05}, ...
    'SelectionFcn', @selectiontournament, ...
    'EliteCount', 2, ...
    'ConstraintTolerance', 1e-6, ...
    'Display', 'iter', ...
    'PlotFcn', {@gaplotbestf, @gaplotmaxconstr}, ...
    'UseParallel', false);

% Run GA
[solution, objectiveValue] = ga(objfun, nItems, [], [], [], [], lb, ub, confun, IntCon, options);

% Display results
disp('Solution:');
disp(solution);
disp(['Total Value: RM ', num2str(-objectiveValue)]);
disp(['Total Weight: ', num2str(sum(weights .* solution))]);
disp(['Total Volume: ', num2str(sum(volumes .* solution))]);

function [c, ceq] = knapsackConstraints(x, weights, volumes, maxWeight, maxVolume, maxItems)
    % Ensure x is a binary vector
    x = round(x);
    c = [sum(weights .* x) - maxWeight; ...
         sum(volumes .* x) - maxVolume; ...
         sum(x) - maxItems];
    ceq = [];
end
