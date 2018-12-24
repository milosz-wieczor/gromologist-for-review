class Entry:
    def __init__(self, content, subsection):
        self.subsection = subsection
        semicol_index = content.find(';')
        if semicol_index >= 0:
            self.content = content[:semicol_index].strip().split()
            self.comment = content[semicol_index:]
        else:
            self.content = content.strip().split()
            self.comment = ''
    
    
class EntryBonded(Entry):
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        self.atoms_per_entry = self.subsection.atoms_per_entry
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
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        self.atoms_per_entry = self.subsection.atoms_per_entry
        self.types = tuple(self.content[:self.atoms_per_entry])
        self.params = None
        self.interaction_type = content[self.atoms_per_entry]
        
        
class EntryAtom(Entry):
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        self.num, self.type, self.resid, self.resname, self.atomname, _, self.charge, self.mass = self.content[:8]
        self.type_b, self.charge_b, self.mass_b = None, None, None
