from itertools import product
from functools import reduce
from .Subsection import *


class Section:
    """
    "Section" is intended to hold e.g. an entire molecule,
    a full set of FF parameters etc.; it should wrap several
    Subsections together
    """
    
    def __init__(self, content, top):
        self.name = 'System'
        self.top = top
        self.dih_processed = False
        self.subsections = [self.yield_sub(content) for content in self.split_content(content)]
    
    def __repr__(self):
        return "{} section with {} subsections".format(self.name, len(self.subsections))
    
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
        # TODO see also Subsection constructor
        special_lines = [n for n, l in enumerate(content) if l.strip().startswith('[')] + [len(content)]
        return [content[beg:end] for beg, end in zip(special_lines[:-1], special_lines[1:])]
        
    def yield_sub(self, content):
        """
        A wrapper that will select which kind of subsection
        should be instantiated (generic, bonded, or params);
        the [ dihedrals ] section gets special treatment as
        first occurrence contains 'proper' and second contains
        'improper' dihedrals, hence we replace to avoid confusion
        :param content: list of strings, content of the subsection
        :return: a Subsection instance (or a derived class)
        """
        until = content[0].index(']')
        header = content[0][:until].strip().strip('[]').strip()
        if header == 'dihedrals':
            if not self.dih_processed:
                self.dih_processed = True
                return SubsectionBonded(content, self)
            else:
                return SubsectionBonded([l.replace('dihedrals', 'impropers') for l in content], self)
        elif header in {'bonds', 'pairs', 'angles', 'settles', 'exclusions', 'cmap', 'position_restraints'}:
            return SubsectionBonded(content, self)
        elif header == 'atoms':
            return SubsectionAtom(content, self)
        elif header == 'moleculetype':
            return SubsectionHeader(content, self)
        elif header in {'defaults', 'atomtypes', 'pairtypes', 'bondtypes', 'angletypes', 'dihedraltypes',
                        'implicit_genborn_params', 'cmaptypes', 'nonbond_params', 'constrainttypes'}:
            return SubsectionParam(content, self)
        else:
            return Subsection(content, self)
        
    def get_subsection(self, section_name):
        """
        Returns the specified subsection; we always need to run merge()
        on SectionParam first to avoid duplicates
        # TODO need special treatment for param sections with different interaction types (mostly dihedraltypes)
        :param section_name:
        :return:
        """
        ssect = [s for s in self.subsections if s.header == section_name]
        if len(ssect) == 0:
            raise KeyError
        elif len(ssect) > 1:
            raise RuntimeError("Error: subsection {} duplicated in {}".format(section_name, str(self)))
        return ssect[0]


class SectionMol(Section):
    """
    This class should wrap the subsections of a single molecule
    (i.e. one [ moleculetype ], one [ atoms ], one [ bonds ] etc.)
    """
    
    def __init__(self, content_list, top):
        self.natoms = None
        self.charge = None
        super().__init__(content_list, top)
        self.bonds = None
        self.mol_name = self.get_subsection('moleculetype').molname
        self.name = '{} molecule'.format(self.mol_name)
    
    def offset_numbering(self, offset, startfrom=0):
        """
        Offsets atom numbering starting from a specified position;
        necessary e.g. when adding or removing atoms to the topology
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        offset = int(offset)
        self._offset_atoms(offset, startfrom)
        self._offset_params(offset, startfrom)

    def _offset_atoms(self, offset, startfrom):
        """
        Offsets atoms in the [ atoms ] section
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        subsection = self.get_subsection('atoms')
        for entry_num, entry in enumerate(subsection):
            if isinstance(entry, EntryAtom) and entry.num >= startfrom:
                entry.num += offset

    def _offset_params(self, offset, startfrom):
        """
        Offsets atomic numbering in all parameter sections,
        e.g., [ bonds ]
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        for sub_name in [s.header for s in self.subsections if s.header != 'atoms']:
            subsection = self.get_subsection(sub_name)
            for entry_num, entry in enumerate(subsection):
                if isinstance(entry, EntryBonded):
                    entry.atom_numbers = tuple(n + (offset * (n >= startfrom)) for n in entry.atom_numbers)
    
    def _get_bonds(self):
        """
        When explicitly asked to, creates a list of bonds stored as
        ordered tuples of atom numbers
        :return: None
        """
        subsection = self.get_subsection('bonds')
        bond_list = []
        for entry in subsection:
            if isinstance(entry, EntryBonded):
                bond_list.append(entry.atom_numbers)
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
        self._make_bond(anchor_own, anchor_other, other)
        if other is not self:
            self._merge_fields(other)
            self.top.sections.remove(other)
            # the stuff below works but is terribly ugly, we need to have API for manipulating content of Top.system
            system_setup = self.top.sections[-1].get_subsection('molecules')
            system_setup._entries = [e for e in system_setup if other.mol_name not in e]

    def _merge_fields(self, other):
        for subs in ['atoms', 'bonds', 'angles', 'pairs', 'dihedrals', 'impropers', 'cmap']:
            # TODO think of other fields as well? like, position_restraints
            # TODO merge all subsections
            try:
                subsection_other = other.get_subsection(subs)
                subsection_own = self.get_subsection(subs)
                subsection_own.add_entries([str(entry) for entry in subsection_other if entry])
            except KeyError:
                pass
    
    def _make_bond(self, atom_own, atom_other, other):
        self._get_bonds()
        other._get_bonds()
        new_bond = [tuple(sorted([int(atom_own), int(atom_other)]))]
        new_angles = self._generate_angles(other, atom_own, atom_other)
        new_pairs, new_dihedrals = self._generate_14(other, atom_own, atom_other)
        for sub, entries in zip(['bonds', 'pairs', 'angles', 'dihedrals'],
                                [new_bond, new_pairs, new_angles, new_dihedrals]):
            subsection = self.get_subsection(sub)
            subsection.add_entries([EntryBonded(subsection.fstring.format(*entry, subsection.prmtype), subsection)
                                    for entry in entries])

    def _generate_angles(self, other, atom_own, atom_other):
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

    def _generate_14(self, other, atom_own, atom_other):
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
        self.name = 'Parameters'
        self.merge()
    
    def merge(self):
        """
        If multiple sections (e.g. [ bondtypes ]) are present in the topology,
        this fn merges them into single sections to avoid searching in all instances
        # TODO do not merge dihedrals if have different interaction types?
        :return: None
        """
        subsection_labels = [sub.label for sub in self.subsections]
        duplicated_subsections = list({label for label in subsection_labels if subsection_labels.count(label) > 1})
        for sub in duplicated_subsections:
            subsections_to_merge = [s for s in self.subsections if s.label == sub]
            merged_subsection = reduce(lambda x, y: x+y, subsections_to_merge)
            position = self.subsections.index(subsections_to_merge[0])
            self.subsections.insert(position, merged_subsection)
            for old in subsections_to_merge:
                self.subsections.remove(old)
