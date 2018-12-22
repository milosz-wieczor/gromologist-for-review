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
        self.sect = section
        self.header = content[0].strip().strip('[]').strip() if '[' in content[0] else 'header'
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
    
    def add_entry(self, new_entry, position=-1):
        if position:
            self.entries.insert(position, new_entry)
        else:
            self.entries.append(new_entry)
    
    def set_entry(self, line_number, new_line):
        self.entries[line_number] = new_line
    
    def get_entry(self, line_number):
        return self.entries[line_number]
        
        
class SubsectionBonded(Subsection):
    """
    SubsectionBonded contains a subsection with entries corresponding to bonded terms,
    e.g., bonds or dihedrals; should be included in SectionMol
    """
    n_atoms = {'bonds': 2, 'pairs': 2, 'angles': 3, 'dihedrals': 4, 'cmap': 5, 'settles': 2, 'exclusions': 3}
    
    def __init__(self, content, section):
        super().__init__(content, section)
        self.atoms_per_entry = SubsectionBonded.n_atoms[self.header]
    
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


class SubsectionParam(Subsection):
    def __init__(self, content, section):
        super().__init__(content, section)
