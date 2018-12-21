class Subsection:
    counter = {}
    
    def __init__(self, content, topology):
        """
        Here we want to have:
          - a unique representation of the section (header/ID)
          - a list of entries (TODO: decide whether to wrap them in another class)
          - a binding to the Top class that holds the Section object

        :param content: list of strings, entire content of the section
        """
        self.top = topology
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


class SubsectionBonded(Subsection):
    n_atoms = {'bonds': 2, 'pairs': 2, 'angles': 3, 'dihedrals': 4, 'cmap': 5}
    
    def __init__(self, content, topology):
        super().__init__(content, topology)
        try:
            self.atoms_per_entry = SubsectionBonded.n_atoms[self.header]
        except KeyError:
            self.atoms_per_entry = None
    
    def find_entry(self, atom):
        pass
