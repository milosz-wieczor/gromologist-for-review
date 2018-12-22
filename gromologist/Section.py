from .Subsection import *


class Section:
    """
    "Section" is intended to hold e.g. an entire molecule,
    a full set of FF parameters etc.; it should wrap several
    Subsections together
    """
    
    def __init__(self, content, top):
        self.top = top
        self.subsections = [self.yield_sub(content) for content in self.split_content(content)]
    
    @staticmethod
    def split_content(content):
        """
        Splits a block of text (list of strings passed to the __init__,
        corresponding to the entire content of the given section)
        into a list of blocs, each starting with a [ section_header ]
        :param content: list of strings, content of section
        :return: list of lists of strings, contents of individual subsections
        """
        # TODO second "dihedral" should be renamed as "improper" to avoid confusion
        special_lines = [n for n, l in content if l.strip().startswith('[')] + [len(content) - 1]
        return [content[beg:end] for beg, end in zip(special_lines[:-1], special_lines[1:])]
        
    def yield_sub(self, content):
        """
        A wrapper that will select which kind of subsection
        should be instantiated (generic, bonded, or params)
        :param content: list of strings, content of the subsection
        :return: a Subsection instance (or a derived class)
        """
        header = content[0].strip().strip('[]').strip()
        if header in {'atoms', 'bonds', 'pairs', 'angles', 'dihedrals', 'settles', 'exclusions', 'cmap'}:
            return SubsectionBonded(content, self)
        elif header in {'defaults', 'atomtypes', 'pairtypes', 'bondtypes', 'angletypes', 'dihedraltypes',
                        'implicit_genborn_params', 'cmaptypes', 'nonbond_params'}:
            return SubsectionParam(content, self)
        else:
            return Subsection(content, self)
        
    def get_subsection(self, section_name):
        """
        Returns the specified subsection; we always need to run merge()
        on SectionParam first to avoid duplicates
        # TODO need special treatment for improper dihedrals
        :param section_name:
        :return:
        """
        if section_name not in ["dihedrals", "impropers"]:
            ssect = [s for s in self.subsections if s.header == section_name]
            if len(ssect) == 0:
                raise KeyError
            elif len(ssect) > 1:
                raise RuntimeError("Error: subsection {} duplicated in {}".format(section_name, str(self)))
            return ssect[0]
        else:
            if section_name == "dihedrals":
                return [s for s in self.subsections if s.header == section_name][0]
            else:
                return [s for s in self.subsections if s.header == section_name][1]


class SectionMol(Section):
    """
    This class should wrap the subsections of a single molecule
    (i.e. one [ moleculetype ], one [ atoms ], one [ bonds ] etc.)
    """
    
    def __init__(self, content_list, top):
        super().__init__(content_list, top)
        self.names_to_nums, self.nums_to_types, self.nums_to_names = None, None, None
        
    def get_dicts(self):
        """
        dicts are not always needed and are costly to calculate,
        so only fill in the values when explicitly asked to
        :return: None
        """
        if not all([self.names_to_nums, self.nums_to_types, self.nums_to_names]):
            self.names_to_nums, self.nums_to_types, self.nums_to_names = self.mol_type_nums()
    
    def mol_type_nums(self):
        """
        Provides bindings between atomnumber and atomtype
        and vice versa for each molecule identified in
        the topology
        :return: tuple of dicts, each dict contains molname:(type:num) and
        molname:(num:type) bindings
        """
        names_to_nums = {}
        nums_to_names = {}
        nums_to_types = {}
        atoms = self.get_subsection('atoms')
        for line in atoms:
            lspl = line.split()
            if len(lspl) > 7 and line.lstrip()[0] not in ["'", '[', ';']:
                names_to_nums["{}-{}-{}".format(lspl[3], lspl[2], lspl[4])] = lspl[0]
                nums_to_names[lspl[0]] = "{}-{}-{}".format(lspl[3], lspl[2], lspl[4])
                nums_to_types[lspl[0]] = lspl[1]
                # TODO add optional reading of state B?
        return names_to_nums, nums_to_types, nums_to_names
    
    def offset_numbering(self, offset, startfrom=0):
        """
        Offsets atom numbering starting from a specified position;
        necessary e.g. when adding or removing atoms to the topology
        :param offset:
        :param startfrom:
        :return:
        """
        offset = int(offset)
        self.offset_atoms(offset, startfrom)
        self.offset_params(offset, startfrom)

    def offset_atoms(self, offset, startfrom):
        subsection = self.get_subsection('atoms')
        for linenum, line in enumerate(subsection):
            lspl = line.split()
            if len(lspl) > 7 and not lspl[0].startswith(';') and int(lspl[0]) >= startfrom:
                # TODO make format string a subsection attribute
                new_line = ('{:6d}{:>11s}{:>7s}{:>7s}{:>7s}{:>7d} ' +
                            (len(lspl) - 6) * '{:>10s} ').format(int(lspl[0]) + offset, *lspl[1:5],
                                                                 int(lspl[5]) + offset, *lspl[6:])
                subsection.set_entry(linenum, new_line)

    def offset_params(self, offset, startfrom):
        for sub_name in [s.header for s in self.subsections if s.header != 'atoms']:
            subsection = self.get_subsection(sub_name)
            n_atoms = subsection.atoms_per_entry
            for linenum, line in enumerate(subsection):
                lspl = line.split()
                if len(lspl) > n_atoms and not lspl[0].startswith(';'):
                    new_line = (n_atoms * '{:>5d} ' + '{:>5s} ' +
                                (len(lspl) - n_atoms - 1)*'{:>10s}').format(*[int(x) + (offset * (int(x) >= startfrom))
                                                                              for x in lspl[:n_atoms]], *lspl[n_atoms:])
                    subsection.set_entry(linenum, new_line)


class SectionParam(Section):
    """
    This class should wrap together sections such as [ bondtypes ],
    [ atomtypes ], [ pairtypes ] etc. and have methods designed to
    facilitate the search of matching params
    """
    
    def __init__(self, content_list, top):
        super().__init__(content_list, top)
    
    def merge(self):
        """
        If multiple sections (e.g. [ bondtypes ]) are present in the topology,
        this fn should merge them into single sections and sort entries to avoid
        searching in all instances
        :return:
        """
        # TODO
        pass
