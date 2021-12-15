# Takes a SMILES representation and creates a xyz file

#imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


#Fragment
def get_fragments(mol,name):
    fragment_names = []
    fragments = Chem.GetMolFrags(mol, asMols = True)
    labels = ["A","B","C"]
    for label,fragment in zip(labels,fragments):
        fragment_names.append(name+label)
    
    return fragments, fragment_names

def generate_conformations(fragments, max_confs=20):
    for fragment in fragments:
        rot_bond = rdMolDescriptors.CalcNumRotatableBonds(fragment)
        confs = min(3 + 3*rot_bond,max_confs)
        AllChem.EmbedMultipleConfs(fragment,numConfs=confs)
    
    return fragments

#Make xyz file
def write_xtb_input_file(fragment, fragment_name):
    number_of_atoms = fragment.GetNumAtoms()
    charge = Chem.GetFormalCharge(fragment)
    #charge = -1
    symbols = [a.GetSymbol() for a in fragment.GetAtoms()] 
    for i,conf in enumerate(fragment.GetConformers()):
        file_name = fragment_name+".xyz"
        with open(file_name, "w") as file:
            file.write(str(number_of_atoms)+"\n")
            file.write("title\n")
            for atom,symbol in enumerate(symbols):
                p = conf.GetAtomPosition(atom)
                line = " ".join((symbol,str(p.x),str(p.y),str(p.z),"\n"))
                file.write(line)
            if charge !=0:
                file.write("$set\n")
                file.write("chrg "+str(charge)+"\n")
                file.write("$end")
        break
                
def  write_input_files(mol,name):
    fragments, fragment_names = get_fragments(mol,name)
    fragments = generate_conformations(fragments)
    for fragment, fragment_name in zip(fragments, fragment_names):
       write_xtb_input_file(fragment, fragment_name)


