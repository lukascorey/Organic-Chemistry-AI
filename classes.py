#! usr/bin/env/python3

class molecule_path_node:
    def __init__(self, molec, pathindex, reactions, recent_reac):
        self.molec = molec
        self.pathindex = pathindex
        self.reactions = reactions 
        self.recent_reac = recent_reac
    def __repr__(self):
        return str(self.recent_reac) + " " + str(self.pathindex) + " " + str(self.reactions)

class atom:
    def __init__(self, element, label):
        self.element = element
        if element == "carbon":
            self.bonds = 4
        elif element == "nitrogen": 
            self.bonds = 3
        elif element == "oxygen":
            self.bonds = 2
        elif element == "hydrogen":
            self.bonds = 1
        elif element == "bromine":
            self.bonds = 1
        self.connected_atoms = []
        self.label = label
        pass

    def connect_atom(self, other):
        self.connected_atoms.append(other)
        other.connected_atoms.append(self)
        if len(self.connected_atoms) > 5: 
            print("bad")
        pass

    def disconnect_atom(self, other):
        new_self = []
        new_other = []
        for i in self.connected_atoms:
            if i.label != other.label: 
                new_self.append(i)
        for i in other.connected_atoms: 
            if i.label != self.label:
                new_other.append(i)
        self.connected_atoms = new_self
        other.connected_atoms = new_other
        pass
        

    def __str__(self):
        return self.element + " atom with " + str(self.bonds) + " bonds labeled " + str(self.label) + "\n"
    
    def __repr__(self):
        return self.label

    def __eq__(self, other): 
        if self.element == other.element and self.label == other.label: 
            return True 
        return False
    
    def __ne__(self, other): 
        if self.element != other.element or self.label != other.label: 
            return True
        return False

class molecule: 
    def __init__(self): 
        self.atoms = []
        self.bonds = []
        pass

    def __repr__(self): 
        return "".join(str(self.bonds))

    def add_atom(self, a):
        self.atoms.append(a)
        pass
    def add_atoms(self, a):
        self.atoms = self.atoms + a
        pass

    def remove_atom(self, a):
        new_atoms = []
        for i in self.atoms: 
            if i != a:
                new_atoms.append(i)
        self.atoms = new_atoms
        pass

    def remove_bond(self, bond):
        new_bonds = []
        for i in self.bonds: 
            if i != bond and i != (bond[1], bond[0], bond[2]):
                new_bonds.append(i)
        self.bonds = new_bonds
        for i in self.atoms:
            if i == bond[0]:
                i.disconnect_atom(bond[1])
            elif i == bond[1]:
                i.disconnect_atom(bond[0])
    def add_bonds(self, bonds_lst): 
        for i in bonds_lst: 
            self.add_bond(i)

    def add_bond(self, bond):
        self.bonds.append(bond)
        for i in self.atoms: 
            if i == bond[0]: 
                for _ in range(bond[2]):
                    i.connect_atom(bond[1])
        
        
    def replace_bond(self, old, new):
        for i in self.atoms: 
            if i == old[0]:
                # note, disconnect atom totally disconnects
                i.disconnect_atom(old[1])
            elif i == old[1]:
                i.disconnect_atom(old[0])
        for i in self.atoms: 
            if i == new[0]: 
                for _ in range(new[2]):
                    i.connect_atom(new[1])
            if i == new[1]:
                for _ in range(new[2]):
                    i.connect_atom(new[0])
        new_bonds = []
        for i in self.bonds: 
            if i != old and (i[1], i[0], 1) != old: 
                new_bonds.append(i)
            else: 
                new_bonds.append(new)
        self.bonds = new_bonds
        pass

    def __eq__(self, other):
        for i in self.bonds: 
            if i not in other.bonds: 
                return False
        return True 

    def __ne__(self, other):
        for i in self.bonds: 
            if i not in other.bonds: 
                return True
        return False 
