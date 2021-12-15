# Renumbering ring indices to maintain the molecule's structure
# ring_num: Renumbers ring indices to maintain molecule structure


from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

# Checks if the indexed character is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False

# Renumbers rings
def rnumber(Compd):
    index = 1
    Comp = Compd
    C2 = ''
    for i in Comp:
        j = i
        if i == '(':
            index+=1
        elif i == ')':
            index-=1
        if is_number(i) == True:
            ind_str = str(index)
            j = ind_str
            Comp = Comp.replace(i, ind_str, 1)
        C2+= j    
    return C2
   
