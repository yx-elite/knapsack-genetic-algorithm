clc, clear, close all

% Parameters
num_items = 10;
weight_limit = 50; % Example weight limit
volume_limit = 30; % Example volume limit
max_items_selected = 5; % Maximum number of items to select
population_size = 100;
num_generations = 100;

% Generate random data for items (weight, volume, value)
weights = [4, 6, 1, 1, 8, 3, 8, 6, 6, 5];
volumes = [6, 5, 7, 2, 1, 6, 4, 4, 9, 3];
values = [62, 32, 59, 95, 48, 14, 69, 84, 6, 46];

% Initialization
population = randi([0, 1], population_size, num_items);

for generation = 1:num_generations
    % Evaluation
    fitness = evaluate_population(population, weights, volumes, values, weight_limit, volume_limit, max_items_selected);
    
    % Selection (Tournament selection)
    selected_indices = tournament_selection(fitness, population_size);
    selected_population = population(selected_indices, :);
    
    % Crossover (Single-point crossover)
    crossover_point = randi([1, num_items - 1]);
    offspring = crossover(selected_population, crossover_point);
    
    % Mutation (Bit-flip mutation)
    mutation_rate = 0.1; % Example mutation rate
    mutated_offspring = mutate(offspring, mutation_rate);
    
    % Replacement
    population = mutated_offspring;
end

% Find the best solution
best_fitness = evaluate_population(population, weights, volumes, values, weight_limit, volume_limit, max_items_selected);
[best_fitness, best_index] = max(best_fitness);
best_solution = population(best_index, :);

disp('Best solution:');
disp(best_solution);
disp('Total value:');
disp(best_fitness);

% Function to evaluate population
function fitness = evaluate_population(population, weights, volumes, values, weight_limit, volume_limit, max_items_selected)
    num_individuals = size(population, 1);
    fitness = zeros(num_individuals, 1);
    for i = 1:num_individuals
        total_weight = sum(population(i, :) .* weights);
        total_volume = sum(population(i, :) .* volumes);
        total_value = sum(population(i, :) .* values);
        if total_weight <= weight_limit && total_volume <= volume_limit && sum(population(i, :)) <= max_items_selected
            fitness(i) = total_value;
        else
            fitness(i) = 0; % Penalize infeasible solutions
        end
    end
end

% Function for tournament selection
function selected_indices = tournament_selection(fitness, population_size)
    tournament_size = 5; % Example tournament size
    selected_indices = zeros(population_size, 1);
    for i = 1:population_size
        tournament_indices = randperm(length(fitness), tournament_size);
        [~, winner_index] = max(fitness(tournament_indices));
        selected_indices(i) = tournament_indices(winner_index);
    end
end

% Function for single-point crossover
function offspring = crossover(selected_population, crossover_point)
    num_parents = size(selected_population, 1);
    offspring = zeros(size(selected_population));
    for i = 1:2:num_parents
        parent1 = selected_population(i, :);
        parent2 = selected_population(i+1, :);
        offspring(i, :) = [parent1(1:crossover_point), parent2(crossover_point+1:end)];
        offspring(i+1, :) = [parent2(1:crossover_point), parent1(crossover_point+1:end)];
    end
end

% Function for bit-flip mutation
function mutated_offspring = mutate(offspring, mutation_rate)
    num_offspring = size(offspring, 1);
    num_genes = size(offspring, 2);
    mutated_offspring = offspring;
    for i = 1:num_offspring
        for j = 1:num_genes
            if rand() < mutation_rate
                mutated_offspring(i, j) = ~mutated_offspring(i, j); % Flip the bit
            end
        end
    end
end
