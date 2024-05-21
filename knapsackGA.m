% EPE472 - ARTIFICIAL INTELLIGENCE AND DATA MINING
% GENETIC ALGORITHM GROUP PROJECT

% LAN YI XIAN (151589)
% LEE JIA MIN (151272)
% LOKE YAN KUANG (152255)

clc, clear, close all


%% Define Knapsack Problem Parameters

% Define data for each items
weights = [5, 1, 9, 9, 6, 5, 2, 5, 6, 7];
volumes = [12, 12, 17, 14, 15, 12, 19, 12, 20, 13];
values = [62.30, 32.53, 59.38, 65.12, 48.67, 14.36, 69.52, 84.71, 6.58, 46.99];
nItems = length(weights);

% Define three constraints
maxWeightLimit = 30;
maxVolumeLimit = 60;
maxItemsSelected = 5;

% Define general parameters
populationSize = 100;
maxGeneration = 100;
tournamentSize = 2;
crossoverRate = 0.7;
mutationRate = 0.01; 

% Initialize population 10 bits (0 for unused, 1 for used)
population = randi([0, 1], populationSize, nItems);

% Initialize fitness and distance plot
[avgFitnessArray, bestFitnessArray, avgDistanceArray] = deal(zeros(maxGeneration, 1));
f = figure;
f.Position(4) = 500;

%% Main Loop of Ginetic Algorithm

for generation = 1:maxGeneration
    % Fitness Evaluation
    fitness = knapsackFitness(population, weights, volumes, values, maxWeightLimit, maxVolumeLimit, maxItemsSelected);

    avgFitnessArray(generation) = mean(fitness);
    bestFitnessArray(generation) = max(fitness);
    avgDistanceArray(generation) = mean(abs(diff(fitness)));
    updatePlot(generation, avgFitnessArray, bestFitnessArray, avgDistanceArray);

    % Selection (Tournament selection)
    selectedPopIndex = tournamentSelection(fitness, populationSize, tournamentSize);
    %selectedPopIndex = rouletteWheelSelection(fitness, populationSize);
    selectedPopulation = population(selectedPopIndex, :);
    
    % Crossover (Single-point crossover)
    %offspring = singlePointCrossover(selectedPopulation, nItems, crossoverRate);
    offspring = twoPointCrossover(selectedPopulation, nItems, crossoverRate);
    
    % Mutation (Bit-flip mutation)
    mutatedOffspring = mutation(offspring, mutationRate);
    population = mutatedOffspring;
end

% Final evaluation to achieve best solution
mutatedFitness = knapsackFitness(mutatedOffspring, weights, volumes, values, maxWeightLimit, maxVolumeLimit, maxItemsSelected);
[bestFitnessVal, bestFitnessIndex] = max(mutatedFitness);
best_solution = population(bestFitnessIndex, :);

fprintf('==============================================\n')
fprintf('     KNAPSACK PROBLEM (GENETIC ALGORITHM)     \n')
fprintf('==============================================\n\n')
fprintf('Best solution \t\t: %s\n', mat2str(best_solution));
fprintf('Total weight \t\t: %d kg\n', sum(mutatedOffspring(1, :) .* weights));
fprintf('Total volume \t\t: %d m^3\n\n', sum(mutatedOffspring(1, :) .* volumes));
fprintf('==============================================\n')
fprintf('Maximised Value \t: RM %.2f\n', bestFitnessVal);
fprintf('==============================================\n')


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
            % Penalize infeasible solutions
            fitness(i) = 0;
        end
    end
end

%% Function for Tournament Selection
function selectedPopIndex = tournamentSelection(fitness, populationSize, tournamentSize)
    selectedPopIndex = zeros(populationSize, 1);
    for i = 1:populationSize
        tournamentIndex = randperm(length(fitness), tournamentSize);
        [~, winnerIndex] = max(fitness(tournamentIndex));
        selectedPopIndex(i) = tournamentIndex(winnerIndex);
    end
end

%% Function for Single-point Crossover Operator
function offspring = singlePointCrossover(selectedPopulation, nItems, crossoverRate)
    nParents = size(selectedPopulation, 1);
    offspring = zeros(size(selectedPopulation));
    for i = 1:2:nParents
        if rand < crossoverRate
            crossoverPoint = randi([1, nItems - 1]);
            parent1 = selectedPopulation(i, :);
            parent2 = selectedPopulation(i + 1, :);
            offspring(i, :) = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
            offspring(i + 1, :) = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
        else
            offspring(i, :) = selectedPopulation(i, :);
            offspring(i + 1, :) = selectedPopulation(i + 1, :);
        end
    end
end

%% Function for Two-point Crossover Operator
function offspring = twoPointCrossover(selectedPopulation, nItems, crossoverRate)
    nParents = size(selectedPopulation, 1);
    offspring = zeros(size(selectedPopulation));
    for i = 1:2:nParents
        if rand < crossoverRate
            crossoverPoints = sort(randi([1, nItems - 1], 1, 2));
            parent1 = selectedPopulation(i, :);
            parent2 = selectedPopulation(i + 1, :);
            offspring(i, :) = [parent1(1:crossoverPoints(1)), parent2(crossoverPoints(1)+1:crossoverPoints(2)), parent1(crossoverPoints(2)+1:end)];
            offspring(i + 1, :) = [parent2(1:crossoverPoints(1)), parent1(crossoverPoints(1)+1:crossoverPoints(2)), parent2(crossoverPoints(2)+1:end)];
        else
            offspring(i, :) = selectedPopulation(i, :);
            offspring(i + 1, :) = selectedPopulation(i + 1, :);
        end
    end
end

%% Function for Mutation Operator
function mutatedOffspring = mutation(offspring, mutationRate)
    nOffsprings = size(offspring, 1);
    nGenes = size(offspring, 2);
    mutatedOffspring = offspring;
    for i = 1:nOffsprings
        if rand < mutationRate
            nMutations = randi([1, nGenes]);
            mutationPoints = randperm(nGenes, nMutations);
            for j = 1:nMutations
                mutationPoint = mutationPoints(j);
                mutatedOffspring(i, mutationPoint) = ~mutatedOffspring(i, mutationPoint);
            end
        end
    end
end

%% Function for Plot Settings
function updatePlot(generation, avgFitnessArray, bestFitnessArray, avgDistanceArray)
    x = 1: generation;
    y1 = avgFitnessArray(1:generation);
    y2 = bestFitnessArray(1:generation);
    y3 = avgDistanceArray(1:generation);
    tiledlayout(2, 1);

    ax1 = nexttile;
    plot(ax1, x, y1, 'b.', 'MarkerSize', 10);
    hold on;
    plot(ax1, x, y2, 'k.', 'MarkerSize', 10);
    title(ax1, ['Best: ', num2str(bestFitnessArray(generation)), '  Mean: ', num2str(avgFitnessArray(generation))]);
    xlabel(ax1, 'Generation')
    ylabel(ax1, 'Fitness Value')
    legend(ax1, 'Mean Fitness', 'Best Fitness', 'Location', 'southeast')

    ax2 = nexttile;
    plot(ax2, x, y3, 'b.', 'MarkerSize', 10);
    title(ax2, 'Average Distance between Individuals');
    xlabel(ax2, 'Generation')
    ylabel(ax2, 'Average Distance')
    drawnow;
end

%% Function for Roulette Wheel Selection
function selectedPopIndex = rouletteWheelSelection(fitness, populationSize)
    totalFitness = sum(fitness);
    selectionProbabilities = fitness / totalFitness;

    % Perform roulette wheel selection
    selectedPopIndex = zeros(populationSize, 1);
    for i = 1:populationSize
        cumulativeProbability = 0;
        for j = 1:length(fitness)
            cumulativeProbability = cumulativeProbability + selectionProbabilities(j);
            if rand <= cumulativeProbability
                selectedPopIndex(i) = j;
                break;
            end
        end
    end
end

    
