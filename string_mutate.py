# Mutation function: Mutation operator Takes 1 molecule as SMILES and replaces one atom/functional group with another.

import random

import Multipoint as co

from rdkit import rdBase


rdBase.DisableLog('rdApp.error')

# List of possible mutations:
def get_symbols():
    symbols = ['C(O)', 'C(OC)', 'C(C)', 'C(F)', 'C(SC)', 'C(C=O)', 'C(C(F)(F)F)', 'C(C#N)',
               'C(N(=O)=O)','C(C4=CC=CC=C4)','C','C(C(=O)OC)','C(N(C)C)','C(N(C))','C(N)']
    return symbols

# Mutation operator:
def mutate(child,mutation_rate):
    if random.random() > mutation_rate:
        return child
    symbols = get_symbols()

    child = list(child)

    for i in range(50):
        mutated_gene = random.randint(0, len(child) - 1)
        random_symbol_number = random.randint(0, len(symbols) - 1)
        new_child = list(child)
        random_number = random.random()

        if new_child[mutated_gene] == 'C':

            # -C-C- mutation
            if new_child[mutated_gene + 1] == 'C':
                new_child[mutated_gene] = symbols[random_symbol_number]
                new_child = co.list2string(new_child)

                if co.string_OK(new_child):
                    return new_child

                return new_child

            # -C-C= mutation
            elif new_child[mutated_gene + 1] == '=' and new_child[mutated_gene - 1] == 'C':
                new_child[mutated_gene] = symbols[random_symbol_number]
                new_child = co.list2string(new_child)
                if co.string_OK(new_child):
                    return new_child
                return new_child

            # -C( mutation : replaces (1) with (2)
            elif new_child[mutated_gene + 1] == '(':
                a = len(child)
                count = 0

                for j in range(mutated_gene, a):
                    if child[j] == '(':
                        count = count + 1
                    elif child[j] == ')':
                        count = count - 1
                    if child[j] == ')' and count == 0:
                        break
                new_child = co.list2string(new_child)

                str1 = new_child[mutated_gene: j + 1]
                str2 = symbols[random_symbol_number]

                if len(str1) < 13:
                    new_child = new_child.replace(str1, str2, 1)
                if co.string_OK(new_child):
                    return new_child
                return new_child

            # -C-C# mutation (rings end in numbers e.g. (c1ccccc1)
            elif new_child[mutated_gene + 1].isnumeric():
                if mutated_gene+2 >= len(new_child):
                    new_child.append(symbols[random_symbol_number][1:])
                    new_child = co.list2string(new_child)
                    if co.string_OK(new_child):
                        return new_child
                    return new_child
                elif new_child[mutated_gene + 2] == 'C':
                    new_child[mutated_gene] = symbols[random_symbol_number]
                    new_child = co.list2string(new_child)
                    if co.string_OK(new_child):
                        return new_child
                    return new_child

                # -C-C#= mutation
                elif new_child[mutated_gene + 2] == '=' and new_child[mutated_gene - 1] == 'C':
                    new_child[mutated_gene] = symbols[random_symbol_number]
                    new_child = co.list2string(new_child)
                    if co.string_OK(new_child):
                        return new_child
                    return new_child

                # -C#( mutation : replaces (1) with (2)
                elif new_child[mutated_gene + 2] == '(':
                    a = len(child)
                    count = 0

                    for j in range(mutated_gene, a):
                        if child[j] == '(':
                            count = count + 1
                        elif child[j] == ')':
                            count = count - 1
                        if child[j] == ')' and count == 0:
                            break

                    new_child = co.list2string(new_child)

                    str1 = new_child[mutated_gene: j + 1]
                    str2 = symbols[random_symbol_number]

                    if len(str1) < 13:
                        new_child = new_child.replace(str1, str2, 1)
                    if co.string_OK(new_child):
                        return new_child
                    return new_child
    return co.list2string(child)








