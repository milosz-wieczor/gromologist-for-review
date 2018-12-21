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
        ssect = [s for s in self.subsections if s.header == section_name]
        assert len(ssect) == 1
        return ssect[0]


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
