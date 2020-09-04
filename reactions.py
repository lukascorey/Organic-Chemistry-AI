#! usr/bin/env/python3
from classes import atom, molecule
from copy import deepcopy
from helpers import remove_molecule_duplicates

def can_hydrobrominate(m):
    for i in m.bonds: 
        # if there are two carbons double bonded to each other 
        if i[0].element == "carbon" and i[1].element == "carbon" and i[2] == 2: 
            # if these two carbons arent connected to anything but other carbons 
            # non_carbons1 = False
            # non_carbons2 = False
            # for j in i[0].connected_atoms: 
            #     if j.element == "oxygen" or j.element == "nitrogen":
            #         non_carbons1 = True
            # for j in i[1].connected_atoms:
            #     if j.element == "oxygen" or j.element == "nitrogen":
            #         non_carbons2 = True
            # if (not (non_carbons1 and non_carbons2)) and (len(i[0].connected_atoms) < 4 or len(i[1].connected_atoms) < 4): 
            return True
    return False

def hydrobromination(m):
    #print("start with: ")
    #print(m)
    new_molecules = []
    # for labeling 
    bromines = 0
    for i in m.atoms: 
        if i.element == "bromine": 
            bromines = bromines + 1
    for i in m.bonds: 
        # if there are two carbons double bonded to each other 
        if i[0].element == "carbon" and i[1].element == "carbon" and i[2] == 2: 
            # if these two carbons arent connected to anything but other carbons 
            # non_carbons1 = False
            # non_carbons2 = False
            # for j in i[0].connected_atoms: 
            #     if j.element == "oxygen" or j.element == "nitrogen":
            #         non_carbons1 = True
            # for j in i[1].connected_atoms:
            #     if j.element == "oxygen" or j.element == "nitrogen":
            #         non_carbons2 = True
            # if no non-carbons and not too many carbons, then we can brominate
            #if (not (non_carbons1 and non_carbons2)): 
            new_molecule = deepcopy(m)
            new_molecule.replace_bond((i[0], i[1], 2), (i[0], i[1], 1))
            bromine = atom("bromine", "bromine" + str(bromines))
            new_molecule.add_atom(bromine)
            # check for which is more stable carbo cation that can fit bromine
            if (len(i[0].connected_atoms) > len(i[1].connected_atoms)): #and (not non_carbons1)) or ((not non_carbons1) and (non_carbons2)):
                new_molecule.add_bond((i[0], bromine, 1))
                new_molecules.append(new_molecule)
            elif (len(i[0].connected_atoms) < len(i[1].connected_atoms)): #and (not non_carbons2)) or ((not non_carbons2) and (non_carbons1)):
                new_molecule.add_bond((i[1], bromine, 1))
                new_molecules.append(new_molecule)
            else: 
                bromine2 = atom("bromine", "bromine" + str(bromines))
                new_molecule2 = deepcopy(new_molecule)
                new_molecule.add_bond((i[0], bromine, 1))
                new_molecule2.add_bond((i[1], bromine2, 1))
                new_molecules.append(new_molecule)
                new_molecules.append(new_molecule2)
    changed = 1
    new_molecules = remove_molecule_duplicates(new_molecules)
    returnlist = new_molecules
    # determine further hydrobromination products
    while changed != 0:
        changed = 0
        replace = []
        for i in new_molecules:
            if can_hydrobrominate(i):
                changed = changed + 1
                replace = replace + hydrobromination(i)
            else:
                replace.append(i)
        new_molecules = remove_molecule_duplicates(replace)
        # append new potential products onto returnlist with dupes removed
        returnlist = returnlist + new_molecules
    #print("end with: ")
    #for i in new_molecules: 
    #    print(i)
    return remove_molecule_duplicates(returnlist)

def can_hydrate(m):
    for i in m.bonds: 
        # if there are two carbons double bonded to each other 
        if i[0].element == "carbon" and i[1].element == "carbon" and i[2] == 2: 
            return True
    return False

def hydration(m):
    #print("start with: ")
    #print(m)
    new_molecules = []
    # for labeling 
    oxygens = 0
    hydrogens = 0
    for i in m.atoms: 
        if i.element == "oxygen": 
            if (int)(i.label[-1]) > oxygens: 
                oxygens = (int)(i.label[-1])
    for i in m.atoms: 
        if i.element == "hydrogen":
            if (int)(i.label[-1]) > hydrogens: 
                hydrogens = (int)(i.label[-1])
    for i in m.bonds: 
        # if there are two carbons double bonded to each other 
        if i[0].element == "carbon" and i[1].element == "carbon" and i[2] == 2: 
            new_molecule = deepcopy(m)
            new_molecule.replace_bond((i[0], i[1], 2), (i[0], i[1], 1))
            oxygen = atom("oxygen", "oxygen" + str(oxygens+1))
            hydrogen = atom("hydrogen", "hydrogen" + str(hydrogens+1))
            new_molecule.add_atom(oxygen)
            new_molecule.add_atom(hydrogen)
            # check for which is more stable carbo cation that can fit bromine
            if (len(i[0].connected_atoms) > len(i[1].connected_atoms)):
                new_molecule.add_bond((i[0], oxygen, 1))
                new_molecule.add_bond((oxygen, hydrogen, 1))
                new_molecules.append(new_molecule)
            elif (len(i[0].connected_atoms) < len(i[1].connected_atoms)):
                new_molecule.add_bond((i[1], oxygen, 1))
                new_molecule.add_bond((oxygen, hydrogen, 1))
                new_molecules.append(new_molecule)
            else: 
                oxygen2 = atom("oxygen", "oxygen" + str(oxygens+1))
                hydrogen2 = atom("hydrogen", "hydrogen" + str(hydrogens+1))
                new_molecule2 = deepcopy(new_molecule)
                new_molecule.add_bond((i[0], oxygen, 1))
                new_molecule.add_bond((oxygen, hydrogen, 1))
                new_molecule2.add_bond((i[1], oxygen2, 1))
                new_molecule2.add_bond((oxygen2, hydrogen2, 1))
                new_molecules.append(new_molecule)
                new_molecules.append(new_molecule2)
    changed = 1
    while changed != 0:
        changed = 0
        replace = []
        for i in new_molecules:
            if can_hydrate(i):
                changed = changed + 1
                replace = replace + hydration(i)
            else:
                replace.append(i)
        new_molecules = replace
    #print("end with: ")
    #for i in new_molecules: 
    #    print(i)
    return remove_molecule_duplicates(new_molecules)

def dehydration(m):
    new_molecules = []
    if can_dehydrate(m):
        for i in m.atoms:
            hydrogen_req = False
            carbon_req = False
            if i.element == "oxygen":
                for j in i.connected_atoms: 
                    if j.element == "hydrogen":
                        hydrogen_req = True
                        hydrogen_atom = j
                    if j.element == "carbon":
                        for k in j.connected_atoms: 
                            if k.element == "carbon" and len(k.connected_atoms)<4:
                                carbon_req = True
                                carbon_atom = j
                # if can dehydrate, i is oxygen, hydrogen_atom is h, carbon_atom is j
                if hydrogen_req and carbon_req:
                    #remove H and O
                    m.remove_bond((i, hydrogen_atom, 1))
                    m.remove_atom(hydrogen_atom)
                    m.remove_bond((carbon_atom, i, 1))
                    m.remove_atom(i)
                    # find most substituted
                    maxi = 0
                    for k in carbon_atom.connected_atoms:
                        if k.element == "carbon":
                            if len(k.connected_atoms) > maxi and len(k.connected_atoms) < 4: 
                                maxi = len(k.connected_atoms)
                    done_lst = [] 
                    for k in carbon_atom.connected_atoms:
                        # for each connected atom to the carbon, if its a max substituted carbon do the dehydr
                        if len(k.connected_atoms) == maxi and k.element == "carbon":
                            #prevents duplicates 
                            if k.label not in done_lst:
                                new_molecule = deepcopy(m)
                                for l in m.bonds: 
                                    if (l[1] == k and l[0] == carbon_atom) or (l[0] == k and l[1] == carbon_atom):
                                        replacebond = l
                                new_molecule.replace_bond(replacebond, (replacebond[0], replacebond[1], replacebond[2] + 1))
                                new_molecules.append(new_molecule)
                                done_lst.append(k.label)
                    break
    changed = 1
    while changed != 0:
        changed = 0
        replace = []
        for i in new_molecules:
            if can_dehydrate(i):
                changed = changed + 1
                replace = replace + dehydration(i)
            else:
                replace.append(i)
        new_molecules = replace
    #print("end with: ")
    #for i in new_molecules: 
    #    print(i)
    return remove_molecule_duplicates(new_molecules)
            
def can_dehydrate(m):
    # need OH bonded to carbon bonded to carbon with a hydrogen
    for i in m.atoms: 
        hydrogen_req = False
        carbon_req = False
        if i.element == "oxygen":
            for j in i.connected_atoms: 
                if j.element == "hydrogen":
                    hydrogen_req = True
                if j.element == "carbon":
                    for k in j.connected_atoms: 
                        if k.element == "carbon":
                            carbon_req = True 
        if hydrogen_req and carbon_req:
            return True
    return False
            

