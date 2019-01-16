

class Entry:
    """
    A generic class representing a single line in the topology.
    In an entry, the actual content and the comments are kept
    in two separate variables.
    """
    def __init__(self, content, subsection):
        self.subsection = subsection
        semicol_index = content.find(';')
        if semicol_index >= 0:
            self.content = content[:semicol_index].strip().split()
            self.comment = content[semicol_index:]
        else:
            self.content = content.strip().split()
            self.comment = ''
    
    def __bool__(self):
        if not self.content and not self.comment:
            return False
        return True
    
    
class EntryBonded(Entry):
    """
    This Entry subclass is intended for entries that correspond
    to bonded interaction (bonds, pairs, angles, dihedrals)
    between specific atoms in the topology
    """
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        self.atoms_per_entry = type(self.subsection).n_atoms[self.subsection.header]
        self.atom_numbers = tuple([int(x) for x in self.content[:self.atoms_per_entry]])
        self.interaction_type = content[self.atoms_per_entry]
        # TODO maybe type assignment should only be performed when asked to, i.e. outside of constructor
        self.types_state_a = None
        self.types_state_b = None
        self.params_state_a = None
        self.params_state_b = None
    
    def read_types(self):
        atoms_sub = self.subsection.section.get_subsection('atoms')
        atoms_sub.get_dicts()
        num_to_type_a = atoms_sub.num_to_type
        num_to_type_b = atoms_sub.num_to_type_b
        self.types_state_a = tuple(num_to_type_a[num] for num in self.atom_numbers)
        self.types_state_b = tuple(num_to_type_b[num] for num in self.atom_numbers)
        
        
class EntryParam(Entry):
    """
    This Entry subclass represents a line containing force field
    parameters, e.g. bondtypes, angletypes, cmaptypes, pairtypes etc.
    that map a set of atom types to a set of FF-specific values
    """
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        self.atoms_per_entry = type(self.subsection).n_atoms[self.subsection.header]
        self.types = tuple(self.content[:self.atoms_per_entry])
        self.params = self.content[self.atoms_per_entry + 1:]
        self.interaction_type = self.content[self.atoms_per_entry]
    
    def __repr__(self):
        return "Parameters entry with atomtypes {}, interaction type {} " \
               "and parameters {}".format(self.types,
                                          self.interaction_type,
                                          self.params)
        
        
class EntryAtom(Entry):
    """
    This Entry subclass corresponds to atoms defined in
    the [ atoms ] section of each molecule
    """
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        try:
            self.num, self.type, self.resid, self.resname, self.atomname, _, self.charge, self.mass = self.content[:8]
        except ValueError:
            self.num, self.type, self.resid, self.resname, self.atomname, _, self.charge = self.content[:7]
            self.mass = 0  # TODO get mass as the default from the ffnonbonded 'atomtypes' section
        self.type_b, self.charge_b, self.mass_b = None, None, None
        self.charge, self.mass = float(self.charge), float(self.mass)
