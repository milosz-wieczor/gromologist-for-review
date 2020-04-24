from itertools import product
from functools import reduce

import gromologist as gml


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
        self.subsections = [self._yield_sub(content) for content in self._split_content(content)]
    
    def __repr__(self):
        return "{} section with {} subsections".format(self.name, len(self.subsections))
    
    @staticmethod
    def _split_content(content):
        """
        Splits a block of text (list of strings passed to the __init__,
        corresponding to the entire content of the given section)
        into a list of blocs, each starting with a [ section_header ]
        :param content: list of strings, content of section
        :return: list of lists of strings, contents of individual subsections
        """
        # TODO second "dihedral" should be renamed as "improper" to avoid confusion
        special_lines = [n for n, l in enumerate(content) if l.strip().startswith('[')] + [len(content)]
        return [content[beg:end] for beg, end in zip(special_lines[:-1], special_lines[1:])]
        
    def _yield_sub(self, content):
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
                return gml.SubsectionBonded(content, self)
            else:
                return gml.SubsectionBonded([l.replace('dihedrals', 'impropers') for l in content], self)
        elif header in {'bonds', 'pairs', 'angles', 'settles', 'exclusions', 'cmap', 'position_restraints'}:
            return gml.SubsectionBonded(content, self)
        elif header == 'atoms':
            return gml.SubsectionAtom(content, self)
        elif header == 'moleculetype':
            return gml.SubsectionHeader(content, self)
        elif header in {'defaults', 'atomtypes', 'pairtypes', 'bondtypes', 'angletypes', 'dihedraltypes',
                        'implicit_genborn_params', 'cmaptypes', 'nonbond_params', 'constrainttypes'}:
            return gml.SubsectionParam(content, self)
        else:
            return gml.Subsection(content, self)
        
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
        
    def __repr__(self):
        return self.name
    
    def select_atoms(self, selection_string):
        sel = gml.SelectionParser(self)
        return sel(selection_string)
    
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
            if isinstance(entry, gml.EntryAtom) and entry.num >= startfrom:
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
                if isinstance(entry, gml.EntryBonded):
                    entry.atom_numbers = tuple(n + (offset * (n >= startfrom)) for n in entry.atom_numbers)
    
    def add_atom(self, atom_number, atom_name, atom_type, charge=0.0, resid=None, resname=None, mass=None):
        """
        For convenience, we try to infer as much as possible
        from existing data, so that it is sufficient to pass
        atom number, atom name and atom type to have a working
        example
        :param atom_number: int, new atom index (1-based)
        :param atom_name: str, name of the atom
        :param atom_type: str, type of the atom
        :param charge: float, charge of the atom
        :param resid: int, residue number
        :param resname: str, residue name
        :param mass: float, mass of the atom
        :return: None
        """
        subs_atoms = self.get_subsection('atoms')
        atoms = subs_atoms.entries
        if not resid and not resname:
            if atom_number > 1:
                ref_entry = [e for e in atoms if (isinstance(e, gml.EntryAtom) and e.num == atom_number - 1)][0]
            else:
                ref_entry = [e for e in atoms if isinstance(e, gml.EntryAtom)][0]
            while not resid:
                q = input("By default, atom will be assigned to residue {}{}. Proceed? [y/n]".format(ref_entry.resname,
                                                                                                     ref_entry.resid))
                if q == 'y':
                    resid = ref_entry.resid
                    resname = ref_entry.resname
                elif q == 'n':
                    return
                else:
                    continue
        elif resid and not resname:
            ref_entry = [e for e in atoms if (isinstance(e, gml.EntryAtom) and e.resid == resid)][0]
            resname = ref_entry.resname
        if not mass:
            param_sect = [s for s in self.top.sections if isinstance(s, SectionParam)][0]
            param_entry = [e for e in param_sect.get_subsection('atomtypes').entries if isinstance(e, gml.EntryParam)
                           and e.content[0] == atom_type][0]
            mass = param_entry.content[2]
        fstring = subs_atoms.fstring
        print(fstring.format(atom_number, atom_type, resid, resname, atom_name, atom_number, charge, mass))
        new_entry = gml.EntryAtom(fstring.format(atom_number, atom_type, resid, resname, atom_name, atom_number,
                                             charge, mass), subs_atoms)
        try:
            position = [n for n, a in enumerate(atoms) if isinstance(a, gml.EntryAtom) and a.num == atom_number][0]
        except IndexError:
            last_atom = [a for a in atoms if isinstance(a, gml.EntryAtom)][-1].num
            if atom_number == last_atom + 1:
                atoms.append(new_entry)
            else:
                raise RuntimeError("Last atom number is {}, cannot create atom nr {}".format(last_atom, atom_number))
        else:
            self.offset_numbering(1, atom_number)
            atoms.insert(position, new_entry)
        self.top.recalc_sys_params()
    
    def del_atom(self, atom_number, del_in_pdb=True):
        self._del_atom(atom_number)
        self._del_params(atom_number)
        self.offset_numbering(-1, atom_number)
        self.top.recalc_sys_params()
        if del_in_pdb:
            if self.top.pdb:
                for to_remove in self._match_pdb_to_top(atom_number):
                    self.top.pdb.delete_atom(to_remove)
    
    def _match_pdb_to_top(self, atom_number):
        """
        Returns a list of PDB atom indices (assuming .top matches .pdb)
        that correspond to the specified atom_number in the molecule topology
        :param atom_number: int, atom number in self (1-based)
        :return: list, PDB atom serials (1-based)
        """
        if not self.top.pdb:
            raise ValueError("No PDB object matched to the currently processed topology")
        count = 0
        pdb_atom_indices = []
        for molecule in self.top.system.keys():
            if molecule != self.mol_name:
                count += self.top.get_molecule(molecule).natoms
            else:
                count += atom_number - 1
                pdb_atom_indices.append(self.top.pdb.atoms[count].serial)
        return pdb_atom_indices
        
    def _del_atom(self, atom_number):
        subsect_atoms = self.get_subsection('atoms')
        chosen = [e for e in subsect_atoms.entries if isinstance(e, gml.EntryAtom) and e.num == atom_number][0]
        subsect_atoms.entries.remove(chosen)
    
    def _del_params(self, atom_number):
        for subs in ['bonds', 'angles', 'pairs', 'dihedrals', 'impropers', 'cmap']:
            try:
                subsection = self.get_subsection(subs)
                to_del = []
                for entry in subsection:
                    if isinstance(entry, gml.EntryBonded):
                        if atom_number in entry.atom_numbers:
                            to_del.append(entry)
                for entry in to_del:
                    subsection.entries.remove(entry)
            except KeyError:
                pass
    
    def _get_bonds(self):
        """
        When explicitly asked to, creates a list of bonds stored as
        ordered tuples of atom numbers
        :return: None
        """
        subsection = self.get_subsection('bonds')
        bond_list = []
        for entry in subsection:
            if isinstance(entry, gml.EntryBonded):
                bond_list.append(entry.atom_numbers)
        self.bonds = bond_list
        
    def add_bond(self, first_atom, second_atom):
        """
        This is just an alias for merge_two if bond is intramolecular
        """
        self.merge_two(self, first_atom, second_atom)

    def merge_two(self, other, anchor_own, anchor_other):
        """
        Creates a new bond by either merging two distinct
        molecules (both being part of the same topology)
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
            system_setup.entries = [e for e in system_setup if other.mol_name not in e]
            self.top.read_system_properties()

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
            subsection.add_entries([gml.EntryBonded(subsection.fstring.format(*entry, subsection.prmtype), subsection)
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
    
    def add_ff_params(self, add_section='all'):
        if add_section == 'all':
            subsections_to_add = ['bonds', 'angles', 'dihedrals', 'impropers']
        else:
            subsections_to_add = [add_section]
        for sub in subsections_to_add:
            try:
                subsection = [s for s in self.subsections if s.header == sub][0]
            except IndexError:
                pass
            else:
                subsection.add_ff_params()


class SectionParam(Section):
    """
    This class should wrap together sections such as [ bondtypes ],
    [ atomtypes ], [ pairtypes ] etc. and have methods designed to
    facilitate the search of matching params
    """
    
    def __init__(self, content_list, top):
        super().__init__(content_list, top)
        self.name = 'Parameters'
        self.defines = {}
        self._merge()
        self._get_defines()
    
    def _merge(self):
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
    
    def _get_defines(self):
        for sub in self.subsections:
            for entry in [e for e in sub.entries if not isinstance(e, gml.EntryParam)]:
                if entry.content and entry.content[0] == "#define":
                    self.top.defines[entry.content[1]] = entry.content[2:]