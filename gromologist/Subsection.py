class Subsection:
    counter = {}
    
    def __init__(self, content, section):
        """
        Here we want to have:
          - a unique representation of the section (header/ID)
          - a list of entries (TODO: decide whether to wrap them in another class)
          - a binding to the Top class that holds the Section object

        :param content: list of strings, entire content of the section
        """
        self.section = section
        # TODO header should read 'impropers', write_header should read 'dihedrals'
        self.header = content[0].strip().strip('[]').strip()
        self.write_header = self.header if self.header != 'impropers' else 'dihedrals'
        if self.header in Subsection.counter.keys():
            Subsection.counter[self.header] += 1
        else:
            Subsection.counter[self.header] = 1
        self.id = Subsection.counter[self.header]
        self.entries = [line for line in content if line.strip() and not line.strip().startswith('[')]
    
    def __str__(self):
        """
        As section headers can be repeated, each section is denoted
        by a header and ID (ID corresponding to the consecutive numbering
        of the specific header, e.g. bonds-3 is the third "bonds" section
        :return: str, section label
        """
        return "{}-{}".format(self.header, self.id)
    
    def __len__(self):
        return len(self.entries)
    
    def __iter__(self):
        self.n = 0
        return self
    
    def __next__(self):
        if self.n == len(self):
            raise StopIteration
        n = self.n
        self.n += 1
        return self.entries[n]
    
    def add_entry(self, new_entry, position=None):
        if position:
            position = int(position)
            self.entries.insert(position, new_entry)
        else:
            self.entries.append(new_entry)
    
    def add_entries(self, new_entries_list, position=None):
        if position:
            position = int(position)
            for new_entry in new_entries_list:
                self.entries.insert(position, new_entry)
                position += 1
        else:
            self.entries.extend(new_entries_list)
    
    def set_entry(self, line_number, new_line):
        self.entries[line_number] = new_line
    
    def get_entry(self, line_number):
        return self.entries[line_number]
        
        
class SubsectionBonded(Subsection):
    """
    SubsectionBonded contains a subsection with entries corresponding to bonded terms,
    e.g., bonds or dihedrals; should be included in SectionMol
    """
    n_atoms = {'bonds': 2, 'pairs': 2, 'angles': 3, 'dihedrals': 4, 'impropers': 4,
               'cmap': 5, 'settles': 2, 'exclusions': 3}
    fstrings = {"{:5} " * n_atoms[x] + '\n' for x in ['bonds', 'pairs', 'angles', 'dihedrals', 'impropers']}
    
    def __init__(self, content, section):
        super().__init__(content, section)
        self.atoms_per_entry = SubsectionBonded.n_atoms[self.header]
        self.prmtype = self.check_parm_type()
    
    def sort(self):
        """
        In case we want to sort entries after some are added at the end of the section
        :return: None
        """
        self.entries.sort(key=self.sorting_fn)
    
    def sorting_fn(self, line):
        """
        Comments should go first, then we sort based on first, second,
        ... column of the section
        :param line: str, line to be sorted
        :return: int, ordering number
        """
        if line.strip().startswith(';'):
            return -1
        lspl = [int(x) for x in line.split()[:self.atoms_per_entry]]
        return sum([i * 10**(4*(self.atoms_per_entry - n)) for n, i in enumerate(lspl)])
            
    def add_ff_params_to_entry(self, atom_list, sect_params):
        """
        Given a bonded term (e.g. "21     24     26    5") converts it to atomtypes,
        finds the respective FF parameters and adds them to the bonded entry
        :param atom_list: iterable, contains atom numbers for the entry in question
        :param sect_params: a SectionParam instance that holds all FF params
        :return: None
        """
        # TODO
        pass
    
    def check_parm_type(self):
        """
        Finds number code for interaction type, e.g. CHARMM uses angletype '5' (urey-bradley)
        while Amber uses angletype '1' (simple harmonic)
        :param parmtype: str, name of the parameter
        :return: str, interaction type
        """
        npar = self.atoms_per_entry
        for line in self:
            lspl = line.split()
            if len(lspl) > npar and not line.strip().startswith(';') and not line.strip().startswith('['):
                return lspl[npar]


class SubsectionParam(Subsection):
    def __init__(self, content, section):
        super().__init__(content, section)


class SubsectionAtom(Subsection):
    def __init__(self, content, section):
        super().__init__(content, section)
        self.fstring = "{:6}{:11}{:7}{:7}{:7}{:7}{:11}{:11}   ; " + '\n'
        self.nat = len([e for e in self.entries if len(e.split()) > 6 and not e.strip().startswith(';')])
        self.section.natoms = self.nat
        self.charge = self.section.charge = self.calc_charge()
    
    def calc_charge(self):
        charge = 0
        for line in self.entries:
            lspl = line.split()
            if len(lspl) > 6 and not lspl[0].startswith(';'):
                charge += float(lspl[6])
        return charge