# Phenyl_GA
This code is a genetic algorithm to generate terphenyl derivative molecules which uses a screening function based on LUMO values, and the average charge of aromatic carbons, of the molecules.

Prerequisites:

QChem
Python packages: a. ase b. rdkit
Installing: Qchem can be downloaded from http://www.q-chem.com/. IMPORTANT: Qchem requires a license. All of the above mentioned packages have conda installtions.

Modules:

Multipoint: Crossover operator Takes 2 molecules as SMILES and produces a new molecule with features of both the parent molecules.

string_mutate: Mutation operator Takes 1 molecule as SMILES and replaces one atom/functional group with another.

xyz_gen: Takes a SMILES representation and creates a xyz file

Sc_fn2: Assigns scores to molecules based on their LUMO energies.

string_GA: Generates a population of molecules that evolve in each generation(iteration) with increasing values of LUMO energies.

terphenyl2.smi: File containing terphenyl SMILES, which are used to make the initial population in string_GA.

ring_num: Renumbers ring indices to maintain molecule structure

Running: Run the file string_GA.
