# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 13:03:08 2026
@author: andrearaaschou
Exercise 2
"""
import numpy as np
import functions as f

# Main variables
pop_size = 200
gene_number = 7
max_generations = 100
num_breed = 20
mutation_rate = 0.1
mutation_free = 5

# Create population
pop = f.create_pop(pop_size, gene_number)
pop = f.evaluate_fitness(pop)
pop = f.sort_pop(pop)

# Empty vectors to store results
topfitness = np.zeros(max_generations)
meanfitness = np.zeros(max_generations)

# main loop over all of our functions
for generation in range(max_generations):
    pop = f.reproduce(pop, num_breed)
    pop = f.mutate(pop, mutation_free, mutation_rate)
    pop = f.evaluate_fitness(pop)
    pop = f.sort_pop(pop)
    
    # save statisitcs
    topfitness[generation] = pop["fitness"][0]
    meanfitness[generation] = np.mean(pop["fitness"])

# plot results
f.plot_convergence(topfitness,meanfitness, max_generations)
f.plot_comparison(pop)

print(f"Genes of the best individual: {np.round(pop["genes"][0, :], 2)}")
print(f"Fitness of the best individual: {np.round(pop["fitness"][0], 2)}")

