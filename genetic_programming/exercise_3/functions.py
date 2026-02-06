# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 13:44:08 2026
@author: andrearaaschou
"""
import numpy as np
import  matplotlib.pyplot as plt
import sys

def target_pol_1(x):
    return -7.5*x**3 + 2*x**2 + 5*x + 2

def target_pol_2(x):
    return 7.5*x**3 - 2*x**2 - 5*x + 2

def target_pol_3(x):
    return -2.5*x**3 * 3*x**2 * (-5.5)*x * (- 2.2)


# Note: This function does not preserve precedence order
# not used
def approx_pol_random(x, genes, chr_size):
    num_op = int(np.ceil(chr_size / 2)) # get the number of operators + 1
    y = genes[chr_size - 1] * x**(num_op - 1) # first step of the polynomial 
    print(f"chr_size is {chr_size}")
    # loop through the rest of the polynomial, get coefficients, operator and power for each step
    for op_i, i in zip(range(chr_size - 2, 0, -2), range(num_op - 2, -1, -1)):
        
        
        a_i = genes[op_i - 1]
        op_i = genes[op_i]
        print(f"a_i is {a_i}")
        print(f"op_i is {op_i}")
        if a_i > 11 or a_i < -11 or op_i < -11 or op_i > 11:
            sys.exit(73)
        """
        if op_i <= 0:
            r = np.random.normal(-5, 2.5)
            if op_i > r:
                y = y + a_i * x ** i
            else:
                y = y * a_i * x ** i
        else:
            r = np.random.uniform(5, 2.5)
            if op_i > r:
                y = y - a_i * x ** i
            else:
                y = y + a_i * x ** i
        """
        op_i += np.random.normal(0, 3.0)
        if op_i <= 0:               # operator is *
            y = y * a_i * x**i
        elif op_i > -5 and op_i < 5: # operator is +
            y = y + a_i * x**i
        else:                        # operator is -
            y = y - a_i * x**i
    return y

def approx_pol(x, genes, chr_size):
    num_op = int(np.ceil(chr_size / 2)) # get the number of operators + 1
    y = genes[chr_size - 1] * x**(num_op - 1) # first step of the polynomial 
    print(f"chr_size is {chr_size}")
    # loop through the rest of the polynomial, get coefficients, operator and power for each step
    #print(f"{len(range(chr_size - 2, 0, -2))} {len(range(num_op - 2, -1, -1))}")
    for op_i, i in zip(range(chr_size - 2, 0, -2), range(num_op - 2, -1, -1)):
        print(f"op_i position is {op_i}")
        a_i = genes[op_i - 1]
        op_i = genes[op_i]
        print(f"a_i is {a_i}")
        print(f"op_i is {op_i}")
        print(f"i is {i}")
        if a_i > 11 or a_i < -11 or op_i < -11 or op_i > 11:
            sys.exit(73)
        if op_i <= -5:               # operator is *
            y = y * a_i * x**i
        elif op_i > -5 and op_i < 5: # operator is +
            y = y + a_i * x**i
        else:                        # operator is -
            y = y - a_i * x**i
    
    return y

# Creates a random population of pop_size chromosomes
# change to allow for different sized chromosomes: if gene is set to 100 it is an empty position.
def create_pop(pop_size, gene_number, pos_gene_numbers):
    genes = np.random.uniform(-10, 10, (pop_size, gene_number))
    fitness = np.zeros((pop_size))
    chr_size = np.random.choice(pos_gene_numbers, size = pop_size) # array with 50 numbers randomly picked from pos_gene_numbers

    for ind in range(pop_size):             # for each individual
        genes[ind, chr_size[ind]:] = 100    # fill the empty part of the array with value = 100

    pop = { # create pouplation (dict)
    "size": pop_size,
    "L": gene_number,
    "chr_size": chr_size,
    "genes": genes,
    "fitness": fitness}
    return pop

# Evaluate fitness for each individual in pop
def evaluate_fitness(pop):
    x = np.linspace(-1, 1, 21)    # x-values (-1,1)
    O = target_pol_2(x)           # observed values
    
    for ind in range(pop["size"]):  # for each individual in pop
        Y = approx_pol(x, pop["genes"][ind], pop["chr_size"][ind])  # send both x and the gene for current individual 
        # Use RMSE formula to calculate fitness and  add another objective: lower chr_size gives higher fitness
        pop["fitness"][ind] = np.sqrt(np.mean((O-Y)**2)) + 0.01 * pop["chr_size"][ind]
    return pop
    
# Sort the population according to fitness rating - small to large
def sort_pop(pop):
    idx = np.argsort(pop["fitness"]) # find order of fitness 
    pop["fitness"] = pop["fitness"][idx] # sort fitness, chromosome size and genes according to fitness order
    pop["genes"] = pop["genes"][idx] 
    pop["chr_size"] = pop["chr_size"][idx]
    return pop

# Reproduce function using a crossover method 
def reproduce(pop, num_breed):
    offspring = []          # empty list to store offspring
    offspring_chr_size = [] # empty list to store offspring chromosome size
    used_parents = []       # empty list to store position indexes of parents already used
    n_offspring = 0               # counter to keep track of successfull matches 
    
    for parent1 in range(0, pop["size"]):    
        if parent1 in used_parents: # if parent1 has already been used
            continue # go to the next possible parent1
        
        # If not: Find parent2
        chr_size = pop["chr_size"][parent1] # get chr_size of parent1
        parent2 = find_match(pop, parent1, used_parents)
        #parent2 = parent1 + 1
        
        if parent2 == -1: # valid match can not be found, go to next individual
            continue
        # else: valid match found, perform recombination
        crossover = np.random.randint(1, min(chr_size, pop["chr_size"][parent2]) - 1) # generate random int [1, number of operators)
                    
        # Create 2 offspring by crossover, save in offspring list
        offspring.append(np.concatenate(
            (pop["genes"][parent1, :crossover], 
             pop["genes"][parent2, crossover:])))
        offspring.append(np.concatenate(
            (pop["genes"][parent2, :crossover], 
             pop["genes"][parent1, crossover:])))
        #offspring_chr_size.extend([pop["chr_size"][parent2], chr_size])
        offspring_chr_size.extend([chr_size, chr_size])
              
        
        used_parents.append(parent1) # add parent1 and parent2 to list with used parents
        used_parents.append(parent2)
        n_offspring += 2 # keep track of how many matches has been made
        
        if n_offspring >= num_breed: # if the wanted number of offspring has been made, break loop
            break

    pop["genes"][-n_offspring:, :] = np.vstack(offspring)   # overwrite the last individuals of pop with offspring
    pop["chr_size"][-n_offspring:] = np.array(offspring_chr_size)
    # update chr_size for offspring
    return pop

# used by reproduce function to find a match for parent1 with the same chromosome size
# input: pop, index of parent1, array with indexes of used parents
# return: index of parent2 (-1 if match can not be found)
def find_match(pop, parent1, used_parents):
    match_found = False                 # keep track of if a successfull match has been found 
    chr_size = pop["chr_size"][parent1] # get chr_size of parent1
       
    # Go through chr_size array in pop starting after parent1, return index of first match
    for i in range(parent1 + 1, pop["size"]):                         # for the rest of the population
        if pop["chr_size"][i] == chr_size:                            # if chr_size is matching
            if i in used_parents: # if individual at index i has already been used
                continue                                              # skip and look for the next match
            parent2 = i                                               # else: valid match with same chr_size has been found
            match_found = True
            break                                                     # break loop
    
    if match_found is False:    # if no match can be found:
        return -1               # return -1 to signal that no match is found
    return parent2              # else if match is found: return index of parent2
    

# Mutate the whole popuolation except the individuals with highest fitness
# try to introduce mutation of chromosome length
def mutate(pop, mutation_free, mutation_rate):
    for ind in range(mutation_free, pop["size"]):       # for each individual
        if np.random.uniform() <= mutation_rate and False:    # if random number is below mutation_rate
            current_size = pop["chr_size"][ind]
            # Generate integer perturbation from current chromosome size
            # low limit: chr_size - 3, high limit: L - chr_size
            #print(f"current size: {current_size}")
            #print(f"current size - 3: {current_size - 3}")
            #print(f"21 - current size : {21 - current_size}")
            var = np.random.randint(-3, 3)
            new_size = max(min(current_size + var, pop["L"]), 3)
            
            
            pop["chr_size"][ind] = new_size
            
            if current_size < new_size: # generate new genes for positions that are not filled from before
                pop["genes"][ind, current_size:new_size] = 0.0
            if current_size > new_size: # fill removed genes with 100
                pop["genes"][ind, new_size:] = 100
            
    
        for gene_ind in range(pop["chr_size"][ind]):    # for each gene in every chromosome 
            if np.random.uniform() <= mutation_rate:    # if random number is below mutation_rate
                pop["genes"][ind, gene_ind] += np.random.normal(0, 0.5) # mutate current gene, sd = 1
                
                # check for values that are out of range (-10, 10) 
                if np.isclose(pop["genes"][ind, gene_ind], 100):
                    print("100 is mutated")
                    print(f"ind is {ind} (individual of pop")
                    print(f"gene_index is {gene_ind}")
                    print(f"chr_size of curren tindividual is {pop["chr_size"][ind]}")
                    continue
                if pop["genes"][ind, gene_ind] > 10:
                    pop["genes"][ind, gene_ind] = 10
                elif pop["genes"][ind, gene_ind] < -10:
                    pop["genes"][ind, gene_ind] = -10
    return pop        

# Make plot with comparison between observed values and predicted values (best individual)
def plot_comparison(pop): 
    x = np.linspace(-1, 1, 21)
    O = target_pol_2(x)
    Y = approx_pol(x, pop["genes"][0], pop["chr_size"][0]) 
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
    
def decode_operator(op):
    if op <= -5:
        return "*"
    elif op < 5:
        return "+"
    else:
        return "-"

def gene_to_expression(genes, chr_size, var="x"):
    genes = genes[:chr_size]
    num_terms = (chr_size + 1) // 2  # number of coefficients

    power = num_terms - 1
    expr = f"{genes[0]:.3g}*{var}^{power}"

    gene_idx = 1
    power -= 1

    while gene_idx + 1 < chr_size:
        op_val = genes[gene_idx]
        coef = genes[gene_idx + 1]

        op = decode_operator(op_val)

        if power > 1:
            term = f"{coef:.3g}*{var}^{power}"
        elif power == 1:
            term = f"{coef:.3g}*{var}"
        else:
            term = f"{coef:.3g}"

        expr += f" {op} {term}"

        gene_idx += 2
        power -= 1

    return expr



