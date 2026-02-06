# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 13:02:29 2026
@author: andrearaaschou
"""
import numpy as np
import  matplotlib.pyplot as plt

def target_pol_1(x):
    return -7.5*x**3 + 2*x**2 + 5*x + 2

def target_pol_2(x):
    return 7.5*x**3 - 2*x**2 - 5*x + 2

def target_pol_3(x):
    return -2.5*x**3 * 3*x**2 * (-5.5)*x * (- 2.2)


def approx_pol(x, genes):
    y = genes[0] * x**3 # first step of the polynomial
    
    # loop through the rest of the polynomial, get coefficients, operator and power for each step
    for a_i, op, i in zip(genes[1:4], genes[4:], [2, 1, 0]): 
        if op <= -5:               # operator is *
            y = y * a_i * x**i
        elif op > -5 and op < 5:   # operator is +
            y = y + a_i * x**i
        else:                      # operator is -
            y = y - a_i * x**i
    return y

# Creates a random population of pop_size chromosomes
def create_pop(pop_size, gene_number):
    genes = np.random.uniform(-10, 10,(pop_size, gene_number)) 
    fitness = np.zeros((pop_size))
    
    pop = { # create pouplation (dict)
    "size": pop_size,
    "L": gene_number,
    "genes": genes,
    "fitness": fitness
    }
    return pop

# Evaluate fitness for each individual in pop
def evaluate_fitness(pop):
    x = np.linspace(-1, 1, 21) # x-values (-1,1)
    O = target_pol_2(x) # observed values
    
    for ind in range(pop["size"]): # for each individual in pop
        Y = approx_pol(x, pop["genes"][ind])  # send both x and the gene for current individual 
        pop["fitness"][ind] = np.sqrt(np.mean((O-Y)**2)) # Use RMSE formula to calculate fitness 
    return pop
    
# Sort the population according to fitness rating - small to large
def sort_pop(pop):
    idx = np.argsort(pop["fitness"]) # find order of fitness 
    pop["fitness"] = pop["fitness"][idx] # sort fitness and genes according to fitness order
    pop["genes"] = pop["genes"][idx]
    return pop

# Reproduce function using a crossover method
def reproduce(pop, num_breed):
    offspring = [] # empty list to store offspring
    
    for parent1 in range(0, num_breed, 2):
        parent2 = parent1 + 1
        crossover = np.random.randint(1, pop["L"]) # generate random int [1, gene_number)
        
        # Create 2 offspring by crossover, save in offspring list
        offspring.append(np.concatenate(
            (pop["genes"][parent1, :crossover], 
             pop["genes"][parent2, crossover:])))
        offspring.append(np.concatenate(
            (pop["genes"][parent2, :crossover], 
             pop["genes"][parent1, crossover:])))
        
    offspring = np.vstack(offspring) # reshape list into 2D array
    pop["genes"][-num_breed:, :]  = offspring # overwrite the last individuals of pop with offspring
    return pop

# Mutate the whole popuolation except the individuals with highest fitness
def mutate(pop, mutation_free, mutation_rate):
    for ind in range(mutation_free + 1, pop["size"]):   # for each individual
        for gene_ind in range(pop["L"]):                # for each gene in every chromosome
            if np.random.uniform() <= mutation_rate:    # if random number is below mutation_rate
                pop["genes"][ind, gene_ind] += np.random.normal() # mutate current gene, sd = 1
                
                # check for values that are out of range (-10, 10)
                if pop["genes"][ind, gene_ind] > 10:
                    pop["genes"][ind, gene_ind] = 10
                elif pop["genes"][ind, gene_ind] < -10:
                    pop["genes"][ind, gene_ind] = -10
    return pop        

# Make plot with comparison between observed values and predicted values (best individual)
def plot_comparison(pop):
    x = np.linspace(-1, 1, 21)
    O = target_pol_2(x)
    Y = approx_pol(x, pop["genes"][0]) 
    plt.plot(x, O, 's', label="Observed values", color="blue")
    plt.plot(x, Y, label="Predicted values", color="red")
    plt.title('Final converged polynomial')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

# Make plot with top and mean fitness for each generation
def plot_convergence(topfitness,meanfitness, max_generations):
    plt.plot(range(max_generations), topfitness, label="Top Fitness", color="blue")
    plt.plot(range(max_generations), meanfitness, label="Mean Fitness", color="black") 
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    plt.title("Convergence plot")
    plt.legend()
    plt.show()
    










