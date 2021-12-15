# GA using SMILES finding terphenyl-based molecules maximizing LUMO
# string_GA: Generates a population of molecules that evolve in each generation(iteration) with increasing values of LUMO energies.
# Imports
from rdkit import Chem
from rdkit import rdBase


rdBase.DisableLog('rdApp.error')

import numpy as np
import random

import Multipoint as co
import string_mutate as mu
import Sc_fn2 as sc
import ring_num as nm


# Read file containing initial population list:
def read_file(file_name):
    smiles_list = []
    if file_name == "Result.out":
        with open(file_name, 'r') as file:
            lines = file.read().splitlines()
            last_line = lines[-5]
            last_line = last_line.strip("[]")
            last_line = last_line.split(",")
            for i in range(len(last_line)):
                last_line[i] = last_line[i].strip("' ")
            smiles_list = last_line
    #Read initial population file
    else:
        with open(file_name, 'r') as file:
            for smiles in file:
                smiles = smiles.replace('\n', '')
                smiles_list.append(smiles)

    return smiles_list


# Make initial population:
def make_initial_population(population_size, file_name):
    smiles_list = read_file(file_name)
    population = []
    if file_name == "Result.out":
        population = smiles_list
    else:
        for _ in range(population_size):
            smiles = random.choice(smiles_list)
            while smiles in population:
                if len(smiles_list) < population_size:
                    break
                smiles = random.choice(smiles_list)
            population.append(smiles)
    return population


# Fitness function:
def calculate_normalized_fitness(scores, charge, population, wt, gaps):
    charges = []
    lumo_scores = []
    HLgaps = []
    for score in scores:
        lumo_scores.append(score)
    for c in charge:
        charges.append(c)
    for g in gaps:
        HLgaps.append(g)
    flat_scores = np.array(lumo_scores)
    flat_scores = flat_scores.flatten()
    sum_lumo = sum(np.absolute(flat_scores))
    charge_scores = np.array(charges)
    charge_scores = charge_scores.flatten()
    sum_charges = sum(np.absolute(charge_scores))
    normalized_scores = [lumo / sum_lumo for lumo in flat_scores]
    normalized_charge = [charge / sum_charges for charge in charge_scores]

    normalized_fitness = []
    corrected_scores = []
    for i in range(len(population)):
        fit = (wt * normalized_scores[i]) + ((1 - wt) * normalized_charge[i]) + (HLgaps[i])
        corrected_scores.append(fit)

    # Scores need normalization because sigma p can't be > 1 in making mating pool
    OldMin = min(corrected_scores)
    OldMax = max(corrected_scores)
    fitness = [(((OldValue - OldMin) * (10 - 1)) / (OldMax - OldMin)) for OldValue in corrected_scores]
    sum_fitness = sum(fitness)
    normalized_fitness = [num / sum_fitness for num in fitness]
    return normalized_fitness


# Makes mating pool where candidates with better fitness have better odds of appearing:
def make_mating_pool(population, fitness, mating_pool_size):
    mating_pool = []
    for i in range(mating_pool_size):
        mating_pool.append(np.random.choice(population, p=fitness))  # Sigma p = 1

    return mating_pool


# Makes new molecules through crossover and mutation operators:
def reproduce(mating_pool, population_size, mutation_rate):
    new_population = []
    while len(new_population) < population_size:
        parent_A = random.choice(mating_pool)
        parent_B = random.choice(mating_pool)
        new_child = co.crossover(parent_A, parent_B)
        while new_child == None:
            new_child = co.crossover(parent_A, parent_B)
        new_child = nm.rnumber(new_child)
        if new_child != None:
            new_child = mu.mutate(new_child, mutation_rate)
            if new_child != None:
                new_population.append(new_child)
    return new_population


# Compares ith generation with (i-1)th generation and returns the best 10 molecules:
def sanitize(population, scores, charges, population_size, wt, gaps):
    smiles_list = []
    population_tuples = []

    fitness = calculate_normalized_fitness(scores, charges, population, wt, gaps)

    for score, string, lumo , charge, gap in zip(fitness, population, scores, charges, gaps):
        canonical_smiles = Chem.MolToSmiles(co.string2mol(string))
        # No duplicates
        if canonical_smiles not in smiles_list:
            smiles_list.append(canonical_smiles)
            population_tuples.append((score, string, lumo,charge, gap))

    population_tuples = sorted(population_tuples, key=lambda x: x[0], reverse=True)[:population_size] #sort by fitness
    new_population = [t[1] for t in population_tuples]
    new_scores = [t[0] for t in population_tuples]
    new_lumo = [t[2] for t in population_tuples]
    new_charges = [t[3] for t in population_tuples]
    new_gaps = [t[4] for t in population_tuples]


    # Prints each generation with respective scores in result.out file
    f = open("Result.out", "a")
    print(new_population, file=f)
    print(new_scores, file=f)
    print(new_lumo, file=f)
    print(new_charges, file=f)
    f.close()

    return new_population, new_lumo, new_charges, new_gaps


# Genetic Algorithm:
def GA(args):
    population_size, file_name, generations, mating_pool_size, mutation_rate, wt = args

    population = make_initial_population(population_size, file_name)

    # Start from initial population
    if file_name != "Result.out":
        f = open("Result.out", "a")
        # Print initial generation:
        print(population, file=f)
        f.close()
        start_generation = 0
        scores, charges, gaps = sc.calculate_scores(population, initial)
    # Restart stopped GA from result file
    else:
        population_size, start_generation, mating_pool_size,charges,scores= getArgs(file_name)


    fitness = calculate_normalized_fitness(scores, charges, population, wt, gaps)
    for generation in range(start_generation,generations):
        mating_pool = make_mating_pool(population, fitness, mating_pool_size)
        new_population = reproduce(mating_pool, population_size, mutation_rate)

        # Print each generation
        f = open("Result.out", "a")
        print(new_population, file=f)
        f.close()

        new_scores, new_charges, new_gaps = sc.calculate_scores(new_population, generation)
        population, scores, charges, gaps = sanitize(population + new_population, scores + new_scores, charges + new_charges,
                                               population_size, wt, gaps + new_gaps)
        fitness = calculate_normalized_fitness(scores, charges, population, wt, gaps)

    return (scores, charges, population, generation + 1, fitness)


# Gets GA args from Result.out to restart previous run.
def getArgs(file_name):
    with open(file_name, 'r') as file:
        lines = file.read().splitlines()
        first_line = lines[0]
        first_line = first_line.strip("[]")
        first_line = first_line.split(",")
        for i in range(len(first_line)):
            first_line[i] = first_line[i].strip("' ")
        initial_smiles = first_line
        count = -1  # To adjust for initial population "50"
        for j in range(len(lines)):
            if lines[j].isnumeric():
                count = count + 1
        charges = lines[-2]
        charges = charges.strip("[]")
        charges = charges.split(",")
        for i in range(len(charges)):
            charges[i] = charges[i].strip(" ")
            charges[i] = float(charges[i])
        LUMO = lines[-3]
        LUMO = LUMO.strip("[]")
        LUMO = LUMO.split(",")
        for i in range(len(LUMO)):
            LUMO[i] = LUMO[i].strip("[] ")
            LUMO[i] = float(LUMO[i])
    generation = count
    population_size = len(initial_smiles)
    next_gen = read_file(file_name)
    mating_pool_size = len(next_gen)

    return (population_size, generation, mating_pool_size,charges, LUMO)


if __name__ == "__main__":
    # Enter GA parameters here

    # 'Result.out' to continue stopped run, terphenyl2.smi to start from terphenyl initial population
    file_name = 'terphenyl2.smi'
    if file_name == "Result.out":
        #get paramaters from results file
        population_size, generation, mating_pool_size,charges,LUMO = getArgs(file_name)
    else:
        population_size = 10
        mating_pool_size = 10
    generations = 10
    mutation_rate = 0.1  # rate at which mutations carried out
    wt = 0.5 # 1.0 weighs only LUMO
    initial = 50  # Just an identifier for initial population

    (scores, charges, population, generation, fitness) = GA([population_size, file_name, generations,
                                                             mating_pool_size, mutation_rate, wt])

    print('done')

    Sl = 0
    for i, j, k in zip(population, scores, charges):
        Sl = Sl + 1
        print(Sl, " Molecule: ", i, "\nScore:", j, "\nCharge:", k)
