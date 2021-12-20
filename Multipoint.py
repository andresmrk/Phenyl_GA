# Crossover function: Crossover operator Takes 2 molecules as SMILES and produces a new molecule with features of both the parent molecules.

from rdkit import Chem

import random

from rdkit import rdBase


rdBase.DisableLog('rdApp.error')

# Checks if the molecule is legitimate:
def string_OK(string):
    mol = string2mol(string)
    #print(mol)
    if not mol:
        return False
    try:
        Chem.SanitizeMol(mol)
        test_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        if test_mol == None:
            return None
        target_size = size_stdev * np.random.randn() + average_size  # parameters set in GA_mol
        print(mol.GetNumAtoms + " -- " + target_size)
        if mol.GetNumAtoms() > 5 and mol.GetNumAtoms() < target_size:
            return True
        else:
            return False
    except:
        return False

def list2string(list):
    string = ''.join(list)
    return string


def string2mol(string):
    try:
        mol = Chem.MolFromSmiles(string)
        return mol
    except:
        return None


# Selects a random splicing point:
def cut_point(parent):
    m = random.randint(0, len(parent) - 1)
    return m


# Multipoint crossover operator:
def crossover(parent_a, parent_b):

    for _ in range(50):
        str1 = ''
        str2 = ''
        # Splicing parent b (gene to be crossed)

        mol = string2mol(str2)
        while mol is None or len(mol.GetAromaticAtoms()) < 6:
            index_b = cut_point(parent_b)
            b = len(parent_b)
            count = 0
            j = index_b
            # Splicing parent_b
            if parent_b[index_b] != '(':
                for j in range(index_b, b):
                    if parent_b[j] == '(':
                        # index_b=1
                        break
            index_b = j

            for j in range(index_b, b):
                if parent_b[j] == '(':
                    count = count + 1
                elif parent_b[j] == ')':
                    count = count - 1
                if parent_b[j] == ')' and count == 0:
                    break
                j = j + 1

            str2 = parent_b[index_b: j + 1]
            mol = string2mol(str2.strip('()'))

        #Crossover with 50/50 chance retention of Parent A genes
        mutated_gene = random.randint(0, len(parent_a)-1)
        new_child = list(parent_a)
        while new_child[mutated_gene] != 'C':
            mutated_gene = random.randint(0, len(parent_a) - 1)
        cross_index = random.random()

        # -C-C- mutation
        if cross_index > 0.5:
            loop = 0
            # swap genes (aromatic rings) by finding a C( or C1(
            while new_child[mutated_gene + 1] != '(':
                mutated_gene = mutated_gene + 1
                # if you reach the end of the molecule, loop from beginning
                if mutated_gene == len(parent_a) - 1:
                    mutated_gene = 0
                    loop = loop + 1
                # no ring groups to swap, add gene to carbon instead
                if loop == 2:
                    mutated_gene = random.randint(0, len(parent_a) - 1)
                    while new_child[mutated_gene] != 'C':
                        mutated_gene = random.randint(0, len(parent_a) - 1)
                    break

        # -C-C- addition
        if new_child[mutated_gene + 1] == 'C':
            new_child[mutated_gene] = 'C' + str2
            child_string = list2string(new_child)
            return child_string

        # -C-C= addition
        elif new_child[mutated_gene + 1] == '=' and new_child[mutated_gene - 1] == 'C':
            new_child[mutated_gene] = 'C' + str2
            child_string = list2string(new_child)
            return child_string

        # -C( cross over rings : replaces (1) with (2)
        elif new_child[mutated_gene + 1] == '(':
            a = len(parent_a)
            count = 0

            for j in range(mutated_gene, a):
                if parent_a[j] == '(':
                    count = count + 1
                elif parent_a[j] == ')':
                    count = count - 1
                if parent_a[j] == ')' and count == 0:
                    break
                # j = j+1

            child_string = list2string(new_child)

            # -C1( mutation
            str1 = child_string[mutated_gene + 1: j + 1]
            mol1 = string2mol(str1.strip('()'))
            if mol1 is not None and len(mol1.GetAromaticAtoms()) >= 6:
                child_string = child_string.replace(str1, str2, 1)
                return child_string
    return child_string

