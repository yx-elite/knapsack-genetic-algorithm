% Main script logic

% Define data for each item
weights = [5, 1, 9, 9, 6, 5, 2, 5, 6, 7];
volumes = [12, 12, 17, 14, 15, 12, 19, 12, 20, 13];
values = [62.30, 32.53, 59.38, 65.12, 48.67, 14.36, 69.52, 84.71, 6.58, 46.99];
nItems = length(weights);

% Define constraints
maxWeightLimit = 40;
maxVolumeLimit = 60;
maxItemsSelected = 5;

% Objective function (negate values to maximize profit)
objectiveFunction = @(x) -sum(values .* x);

% Set the options for the genetic algorithm
options = optimoptions('ga', 'Display', 'iter', 'UseParallel', true);

% Define lower and upper bounds (0 or 1 for each item)
lb = zeros(1, nItems);
ub = ones(1, nItems);

% Solve the problem using ga
[x, fval] = ga(objectiveFunction, nItems, [], [], [], [], lb, ub, @combinedConstraints, options);

% Display results
selectedItems = find(x);
disp('Selected item indices:');
disp(selectedItems);
disp('Total value:');
disp(-fval);
disp('Total weight:');
disp(sum(weights .* x));
disp('Total volume:');
disp(sum(volumes .* x));

% Define combined constraints function
function [c, ceq] = combinedConstraints(x)
    weights = [5, 1, 9, 9, 6, 5, 2, 5, 6, 7];
    volumes = [12, 12, 17, 14, 15, 12, 19, 12, 20, 13];
    maxWeightLimit = 40;
    maxVolumeLimit = 60;
    maxItemsSelected = 5;
    weightConstraint = sum(weights .* x) - maxWeightLimit;
    volumeConstraint = sum(volumes .* x) - maxVolumeLimit;
    itemsConstraint = sum(x) - maxItemsSelected;
    c = [weightConstraint, volumeConstraint, itemsConstraint];
    ceq = [];
end
