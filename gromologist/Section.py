from .Subsection import Subsection


class Section:
    """
    "Section" is intended to hold e.g. an entire molecule,
    a full set of FF parameters etc.; it should wrap several
    Subsections together
    """
    
    def __init__(self):
        self.subsections = []


class SectionMol(Section):
    """
    This class should wrap the subsections of a single molecule
    (i.e. one [ moleculetype ], one [ atoms ], one [ bonds ] etc.)
    """
    
    def __init__(self):
        super().__init__()
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
    
    def get_subsection(self, section_name):
        """
        Returns the specified subsection
        # TODO need special treatment for improper dihedrals
        :param section_name:
        :return:
        """
        ssect = [s for s in self.subsections if s.header == section_name]
        assert len(ssect) == 1
        return ssect[0]


class SectionParam(Section):
    """
    This class should wrap together sections such as [ bondtypes ],
    [ atomtypes ], [ pairtypes ] etc. and have methods designed to
    facilitate the search of matching params
    """
    
    def __init__(self):
        super().__init__()
