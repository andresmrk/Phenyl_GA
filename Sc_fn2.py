# Scoring function (Evaluation function): Sc_fn2: Assigns scores to molecules based on their LUMO energies.
from rdkit import Chem

from ase.calculators.qchem import QChem
from ase.optimize import LBFGS
from ase.io import read

import xyz_gen as xg


# Opens generated Qchem output file and extracts LUMO (followed by equivalent HOMO function):
def extract_LUMO(file):
    number = 0
    with open(file, "r") as output_file:
        for line in output_file:
            if '-- Virtual --' in line:
                number += 1
    # print(number)
    with open(file, "r") as output_file:
        LUMO = []
        count = 0
        for line in output_file:
            if count == number:
                LUMO.append(line.split()[0])
                break
            if '-- Virtual --' in line:
                count += 1
        LUMO = LUMO[0:6]
        for i in range(len(LUMO)):
            LUMO[i] = float(LUMO[i])
    return LUMO[0]

def extract_HOMO(file):
    number = 0
    with open(file, "r") as output_file:
        for line in output_file:
            if '-- Virtual --' in line:
                number += 1
    # print(number)
    with open(file, "r") as output_file:
        HOMO = []
        temp = []
        count = 0
        prevLine = ""
        secondprevLine = ""
        for line in output_file:
            if count == number:
                temp = secondprevLine.split()
                HOMO.append(temp[len(temp)-1])
                break
            if '-- Virtual --' in line:
                count += 1
            secondprevLine = prevLine
            prevLine = line
        for i in range(len(HOMO)):
            HOMO[i] = float(HOMO[i])

    return HOMO[0]

# extracts mulliken charges of all atoms from Qchem output
def extract_charges(file):
    number = 0
    with open(file, "r") as output_file:
        for line in output_file:
            if 'Charge (a.u.)' in line:
                number += 1
    print(number)
    with open(file, "r") as output_file:
        charges = []
        count = 0
        for line in output_file:
            if count == number:
                while "Sum of atomic charges" not in line:
                    line = output_file.readline()
                    # print(line)
                    if "----------------------------------------" in line:
                        break
                    atom_charge = line.split()
                    # print(atom_charge)
                    charges.append(float(atom_charge[2]))
                break
            if 'Charge (a.u.)' in line:
                count += 1
    return charges

def calcHOMOLUMOGap(homo, lumo):
    gap = lumo - homo
    gap = gap * 27.2113961 #Converting from hartree (a.u) to eV
    h = 4.1356677 * 10 ** -15  # eV s (Plancks constant)
    C = 3 * 10 ** 17  # nm /s (Speed of Light)
    wavelength = h * C / gap # from E = hv and C = y*v
    return wavelength

# takes a SMILES string, and gets indeces of all aromatic carbon atoms.
def get_aromatic_indeces(SMILES_string):
    mol = Chem.MolFromSmiles(SMILES_string)
    c_aromatic = []
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetSmarts() == 'c':
            c_aromatic.append(i)
    return c_aromatic


# Function to calculate the scores of the molecules in the population:
def calculate_scores(population, generation):
    scores = []
    charges = []
    gaps = []
    # smiles to xyz
    for i in range(len(population)):
        name = str(generation) + "mol" + str(i)
        mol = Chem.MolFromSmiles(population[i])
        mol = Chem.AddHs(mol)
        xg.write_input_files(mol, name)

        comp = read(name + 'A.xyz')
        calc = QChem(jobtype='opt',
                     label=name,
                     method='PBE',
                     basis='3-21G',
                     scf_convergence='8',
                     scf_max_cycles='500')
        comp.set_calculator(calc)
        opt = LBFGS(comp)
        opt.run()

        # LUMO from qchem.out
        file = name + ".out"
        score = extract_LUMO(file)
        scores.append(score)
        HOMO = extract_HOMO(file)

        #Penalized for being outside 200-400nm
        gap = calcHOMOLUMOGap(HOMO,score)
        if gap >= 200 and gap <= 400:
            gaps.append(0)
        else:
            gaps.append(-1)

        # Mulliken Charge average of aromatic atoms from qchem.out
        charge_list = extract_charges(file)
        charge_indeces = get_aromatic_indeces(population[i])
        charge = []

        #extract only charges of aromatic carbons to average
        for j in range(len(charge_indeces)):
            charge.append(charge_list[charge_indeces[j]])
        if len(charge) == 0:
            mean = 0
        else:
            mean = sum(charge) / len(charge)

        charges.append(mean)

    # Printing scores:
    f = open("Result.out", "a")
    print(generation, file=f)
    print(scores, file=f)
    print(charges, file=f)
    f.close()

    return (scores, charges, gaps)
