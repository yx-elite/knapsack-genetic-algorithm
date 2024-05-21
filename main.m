clc, clear, close all

%% Define Knapsack Problem Parameters

% Define data for each items
weights = [4, 6, 1, 1, 8, 3, 8, 6, 6, 5];
volumes = [6, 5, 7, 2, 1, 6, 4, 4, 9, 3];
values = [62, 32, 59, 95, 48, 14, 69, 84, 6, 46];
nItems = length(weights);

% Define three constraints
maxWeightLimit = 50;
maxVolumeLimit = 30;
maxItemsSelected = 5;

% Define general parameters
populationSize = 100;
maxGeneration = 100;

% Initialize population 10 bits (0 for unused, 1 for used)
population = randi([0, 1], populationSize, nItems);

%% Main Loop of Ginetic Algorithm

for generation = 1:maxGeneration
    % Evaluation
    fitness = knapsackFitness(population, weights, volumes, values, maxWeightLimit, maxVolumeLimit, maxItemsSelected);
    
    % Selection (Tournament selection)
    selectedPopIndex = tournamentSelection(fitness, populationSize);
    selectedPopulation = population(selectedPopIndex, :);
    
    % Crossover (Single-point crossover)
    crossoverPoint = randi([1, nItems - 1]);
    offspring = singlePointCrossover(selectedPopulation, crossoverPoint);
    
    % Mutation (Bit-flip mutation)
    mutationRate = 0.01; 
    mutatedOffspring = mutation(offspring, mutationRate);
    population = mutatedOffspring; 
end

% Find the best solution
mutatedFitness = knapsackFitness(population, weights, volumes, values, maxWeightLimit, maxVolumeLimit, maxItemsSelected);
[bestFitnessVal, bestFitnessIndex] = max(mutatedFitness);
best_solution = population(bestFitnessIndex, :);

disp('Best solution:');
disp(best_solution);
disp('Total value:');
disp(bestFitnessVal);

%% Function for Fitness Evaluation
function fitness = knapsackFitness(population, weights, volumes, values, maxWeightLimit, maxVolumeLimit, maxItemsSelected)
    fitness = zeros(size(population, 1), 1);
    for i = 1:size(population, 1)
        totalWeight = sum(population(i, :) .* weights);
        totalVolume = sum(population(i, :) .* volumes);
        totalValue = sum(population(i, :) .* values);
        if totalWeight <= maxWeightLimit && totalVolume <= maxVolumeLimit && sum(population(i, :)) <= maxItemsSelected
            fitness(i) = totalValue;
        else
            fitness(i) = 0; % Penalize infeasible solutions
        end
    end
end

%% Function for Tournament Selection
function selectedPopIndex = tournamentSelection(fitness, populationSize)
    tournamentSize = 5;
    selectedPopIndex = zeros(populationSize, 1);
    for i = 1:populationSize
        tournamentIndex = randperm(length(fitness), tournamentSize);
        [~, winnerIndex] = max(fitness(tournamentIndex));
        selectedPopIndex(i) = tournamentIndex(winnerIndex);
    end
end

%% Function for Single-point Crossover Operator
function offspring = singlePointCrossover(selectedPopulation, crossoverPoint)
    nParents = size(selectedPopulation, 1);
    offspring = zeros(size(selectedPopulation));
    for i = 1:2:nParents
        parent1 = selectedPopulation(i, :);
        parent2 = selectedPopulation(i+1, :);
        offspring(i, :) = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
        offspring(i+1, :) = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
    end
end

%% Function for Mutation Operator
function mutatedOffspring = mutation(offspring, mutationRate)
    nOffsprings = size(offspring, 1);
    nGenes = size(offspring, 2);
    mutatedOffspring = offspring;
    for i = 1:nOffsprings
        for j = 1:nGenes
            if rand < mutationRate
                mutatedOffspring(i, j) = ~mutatedOffspring(i, j);
            end
        end
    end
end
