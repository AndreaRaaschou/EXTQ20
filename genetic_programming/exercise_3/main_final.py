# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 10:50:42 2026
@author: andrearaaschou
Exercise 3
"""
import numpy as np
import functions_final as f

# Main variables
pop_size = 1000
gene_number = 21
pos_gene_numbers = np.linspace(3, gene_number, num = 10).astype(int) # possible number of genes (corresponding to 1th-10th degree polynomial)
max_generations = 100
num_breed = 2 * int(0.08 * pop_size)
mutation_rate = 0.05
mutation_free = int(0.05 * pop_size)

# Create population
pop = f.create_pop(pop_size, gene_number, pos_gene_numbers)
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

print(f"Genes of the best individual: {np.round(pop["genes"][0, :pop["chr_size"][0]], 2)}")
print(f"Fitness of the best individual: {np.round(pop["fitness"][0], 2)}")
print("Polynomial of the best individual")
print(f.gene_to_expression(pop["genes"][0, :], pop["chr_size"][0]))


