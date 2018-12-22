from .Subsection import *
from itertools import product


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
        elif header == 'atoms':
            return SubsectionAtom(content, self)
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
        self.natoms = None
        super().__init__(content_list, top)
        self.names_to_nums, self.nums_to_types, self.nums_to_names = None, None, None
        self.bonds = None
        
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
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        offset = int(offset)
        self.offset_atoms(offset, startfrom)
        self.offset_params(offset, startfrom)

    def offset_atoms(self, offset, startfrom):
        """
        Offsets atoms in the [ atoms ] section
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        subsection = self.get_subsection('atoms')
        for linenum, line in enumerate(subsection):
            lspl = line.split()
            if len(lspl) > 7 and not lspl[0].startswith(';') and int(lspl[0]) >= startfrom:
                new_line = ('{:6d}{:>11s}{:>7s}{:>7s}{:>7s}{:>7d} ' +
                            (len(lspl) - 6) * '{:>10s} ').format(int(lspl[0]) + offset, *lspl[1:5],
                                                                 int(lspl[5]) + offset, *lspl[6:])
                subsection.set_entry(linenum, new_line)

    def offset_params(self, offset, startfrom):
        """
        Offsets atomic numbering in all parameter sections,
        e.g., [ bonds ]
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
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
    
    def get_bonds(self):
        """
        When explicitly asked to, creates a list of bonds stored as
        ordered tuples of atom numbers
        :return: None
        """
        subsection = self.get_subsection('bonds')
        bond_list = []
        for line in subsection:
            if line.split() and line.split()[2] == "1":
                bond_list.append((int(line.split()[0]), int(line.split()[1])))
        self.bonds = bond_list

    def merge_two(self, other, anchor_own, anchor_other):
        """
        Creates a new bond by either merging two distinct molecules
        or adding a new bond within a single molecule
        :param other: an SectionMol instance, the other molecule that participates in the bond (can be self)
        :param anchor_own: int, number of the atom that will form the new bond in self
        :param anchor_other: int, number of the atom that will form the new bond in other (or self, if other is self)
        :return: None
        """
        anchor_other = int(anchor_other)
        anchor_own = int(anchor_own)
        if other is not self:
            other.offset_numbering(self.natoms)
            anchor_other += self.natoms
        self.make_bond(anchor_own, anchor_other, other)
        self.merge_fields(other)

    def merge_fields(self, other):
        for subs in ['atoms', 'bonds', 'angles', 'pairs', 'dihedrals', 'impropers', 'cmap']:
            # TODO go for try/except
            # TODO need to get rid of other somehow after is transferred to self
            # TODO think of other fields as well?
            subsection_other = other.get_subsection(subs)
            subsection_own = self.get_subsection(subs)
            subsection_own.add_entries([line for line in subsection_other if line.strip()])
    
    def make_bond(self, atom_own, atom_other, other):
        # TODO needs to be wrapped in merge_two
        other.get_bonds()
        new_bond = [tuple(sorted([int(atom_own), int(atom_other)]))]
        new_angles = self.generate_angles(other, atom_own, atom_other)
        new_pairs, new_dihedrals = self.generate_14(other, atom_own, atom_other)
        for sub, entries in zip(['bonds', 'pairs', 'angles', 'dihedrals'],
                                [new_bond, new_pairs, new_angles, new_dihedrals]):
            subsection = self.get_subsection(sub)
            subsection.add_entries([subsection.fstring.format(*entry, subsection.prmtype) for entry in entries])

    def generate_angles(self, other, atom_own, atom_other):
        """
        Generates new angles when an additional bond is formed
        :param other: SectionMol instance, the other molecule that participates in the bond (can be self)
        :param atom_own:
        :param atom_other:
        :return:
        """
        neigh_atoms_1 = [[b for b in bond if b != atom_own][0] for bond in self.bonds if atom_own in bond]
        neigh_atoms_2 = [[b for b in bond if b != atom_other][0] for bond in other.bonds if atom_other in bond]
        new_angles = [(at1, atom_own, atom_other) for at1 in neigh_atoms_1]
        new_angles += [(atom_own, atom_other, at2) for at2 in neigh_atoms_2]
        return new_angles

    def generate_14(self, other, atom_own, atom_other):
        """
        Generates new 1-4 interaction (pairs and dihedrals)
        when an additional bond is formed
        :param other:
        :param atom_own:
        :param atom_other:
        :return:
        """
        # atoms directly neighboring with the new bond
        neigh_atoms_1 = [[b for b in bond if b != atom_own][0] for bond in self.bonds if atom_own in bond]
        neigh_atoms_2 = [[b for b in bond if b != atom_other][0] for bond in other.bonds if atom_other in bond]
        # atoms only neighboring with atoms from the above lists
        neigh_atoms_11 = [list(set(bond).difference(set(neigh_atoms_1)))[0] for bond in self.bonds
                          if set(neigh_atoms_1) & set(bond) and atom_own not in bond]
        neigh_atoms_21 = [list(set(bond).difference(set(neigh_atoms_2)))[0] for bond in other.bonds
                          if set(neigh_atoms_2) & set(bond) and atom_other not in bond]
        new_pairs = list(product(neigh_atoms_1, neigh_atoms_2)) + list(product([atom_own], neigh_atoms_21)) + \
            list(product([atom_other], neigh_atoms_11))
        new_dihedrals = [(a, atom_own, atom_other, d) for a, d in list(product(neigh_atoms_1, neigh_atoms_2))]
        new_dihedrals += [(a, b, atom_own, atom_other) for a in neigh_atoms_11 for b in neigh_atoms_1
                          if (a, b) in self.bonds or (b, a) in self.bonds]
        new_dihedrals += [(atom_own, atom_other, c, d) for d in neigh_atoms_21 for c in neigh_atoms_2
                          if (c, d) in self.bonds or (d, c) in self.bonds]
        return new_pairs, new_dihedrals


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
