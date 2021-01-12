from itertools import product, combinations
from functools import reduce
from copy import deepcopy
from glob import glob

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
                return gml.SubsectionBonded([line.replace('dihedrals', 'impropers') for line in content], self)
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

    @property
    def atoms(self):
        sub = self.get_subsection('atoms')
        return [entry for entry in sub if isinstance(entry, gml.EntryAtom)]
    
    def select_atoms(self, selection_string):
        """
        Returns atoms' indices according to the specified selection string
        :param selection_string: str, a VMD-compatible selection
        :return: list, 0-based indices of atoms compatible with the selection
        """
        sel = gml.SelectionParser(self)
        return sel(selection_string)  # TODO enable getting atoms' properties

    def select_atom(self, selection_string):
        """
        Returns atoms' indices according to the specified selection string
        :param selection_string: str, a VMD-compatible selection
        :return: int, 0-based index of atom compatible with the selection
        """
        sel = gml.SelectionParser(self)
        result = sel(selection_string)
        if len(result) > 1:
            raise RuntimeError("Selection {} returned more than one atom: {}".format(selection_string, result))
        elif len(result) < 1:
            raise RuntimeError("Selection {} returned no atoms".format(selection_string, result))
        return result[0]

    def print_molecule(self):
        sub = self.get_subsection('atoms')
        for entry in sub:
            print(str(entry), end='')

    @property
    def is_alchemical(self):
        sect = self.get_subsection('atoms')
        for ent in sect.entries:
            if isinstance(ent, gml.EntryAtom) and ent.type_b is not None:
                return True
        return False
    
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

    def gen_state_b(self, atomname=None, resname=None, resid=None, atomtype=None, new_type=None, new_charge=None,
                    new_mass=None):
        """
        Generates alchemical state B for a subset of atoms,
        with specified types/charges/masses
        :param atomname: str, these atomnames will be selected
        :param resname: str, these residue names will be selected
        :param resid: int, these residue IDs will be selected
        :param atomtype: str, these atomtypes will be selected
        :param new_type: str, new value for atomtype (default is copy from state A)
        :param new_charge: float, new value for charge (default is copy from state A)
        :param new_mass: float, new value for mass (default is copy from state A)
        :return: None
        """
        sub = self.get_subsection('atoms')
        for entries in [entry for entry in sub if isinstance(entry, gml.EntryAtom)]:
            criteria = all([(atomname is None or entries.atomname == atomname),
                            (resname is None or entries.resname == resname),
                            (resid is None or int(entries.resid) == int(resid)),
                            (atomtype is None or entries.type == atomtype)])
            if criteria:
                entries.type_b = new_type if new_type is not None else entries.type
                entries.mass_b = new_mass if new_mass is not None else entries.mass
                entries.charge_b = new_charge if new_charge is not None else entries.charge

    def drop_state_a(self, remove_dummies=False, atomname=None, resname=None, resid=None, atomtype=None):
        """
        Collapses alchemical B states, making state B
        the new non-alchemical default state A
        :param remove_dummies: bool, whether to remove B-state dummies
        :param atomname: str, name of the selected atom(s) for which state B will be dropped
        :param resname: str, name of the selected residue(s) for which state B will be dropped
        :param resid: int, number of the selected residue(s) for which state B will be dropped
        :param atomtype: str, type of the selected atom(s) for which state B will be dropped
        :return: None
        """
        if not remove_dummies:
            print("Warning: dropping state A parameters, but keeping dummies (if exist). To remove all atoms with "
                  "type names starting with D, rerun this fn with 'remove_dummies=True'.")
        if atomname or resname or resid or atomtype:
            selected = set()
            sub = self.get_subsection('atoms')
            for entries in [entry for entry in sub if isinstance(entry, gml.EntryAtom)]:
                criteria = all([(atomname is None or entries.atomname == atomname),
                                (resname is None or entries.resname == resname),
                                (resid is None or int(entries.resid) == int(resid)),
                                (atomtype is None or entries.type == atomtype)])
                if criteria:
                    selected.add(entries.num)
        else:
            selected = list(range(1, self.natoms+1))
        if remove_dummies:
            sub = self.get_subsection('atoms')
            dummies = [entry for entry in sub if isinstance(entry, gml.EntryAtom) and entry.type_b and
                       entry.type_b[0] == "D" and entry.num in selected]
            print(dummies)
            while dummies:
                to_remove = dummies[0]
                self.del_atom(to_remove.num)
                dummies = [entry for entry in sub if isinstance(entry, gml.EntryAtom) and entry.type_b and
                           entry.type_b[0] == "D" and entry.num in selected]
        for sub in self.subsections:
            for entry in sub:
                if (isinstance(entry, gml.EntryAtom) and entry.num in selected) \
                        or (isinstance(entry, gml.EntryBonded) and any([x in selected for x in entry.atom_numbers])):
                    if isinstance(entry, gml.EntryAtom) and entry.type_b is not None:
                        entry.type, entry.mass, entry.charge = entry.type_b, entry.mass_b, entry.charge_b
                        entry.type_b, entry.mass_b, entry.charge_b = 3 * [None]
                    elif isinstance(entry, gml.EntryBonded) and entry.params_state_b:
                        entry.params_state_a = entry.params_state_b
                        entry.params_state_b = []
                    if isinstance(entry, gml.EntryBonded) and entry.types_state_b is not None:
                        entry.types_state_a = entry.types_state_b
                        entry.types_state_b = None

    def swap_states(self, atomname=None, resname=None, resid=None, atomtype=None):
        """
        Swaps alchemical states A and B
        :param atomname: str, name of the selected atom(s) for which state B will be swapped
        :param resname: str, name of the selected residue(s) for which state B will be swapped
        :param resid: int, number of the selected residue(s) for which state B will be swapped
        :param atomtype: str, type of the selected atom(s) for which state B will be swapped
        :return: None
        """
        if atomname or resname or resid or atomtype:
            selected = set()
            sub = self.get_subsection('atoms')
            for entries in [entry for entry in sub if isinstance(entry, gml.EntryAtom)]:
                criteria = all([(atomname is None or entries.atomname == atomname),
                                (resname is None or entries.resname == resname),
                                (resid is None or int(entries.resid) == int(resid)),
                                (atomtype is None or entries.type == atomtype)])
                if criteria:
                    selected.add(entries.num)
        else:
            selected = list(range(1, self.natoms + 1))
        for sub in self.subsections:
            for entry in sub:
                if (isinstance(entry, gml.EntryAtom) and entry.num in selected) \
                        or (isinstance(entry, gml.EntryBonded) and any([x in selected for x in entry.atom_numbers])):
                    if isinstance(entry, gml.EntryAtom) and entry.type_b is not None:
                        (entry.type, entry.mass, entry.charge, entry.type_b, entry.mass_b, entry.charge_b) = \
                            (entry.type_b, entry.mass_b, entry.charge_b, entry.type, entry.mass, entry.charge)
                    elif isinstance(entry, gml.EntryBonded) and entry.params_state_b:
                        entry.params_state_a, entry.params_state_b = entry.params_state_b, entry.params_state_a
                    if isinstance(entry, gml.EntryBonded) and entry.types_state_b is not None:
                        entry.types_state_a, entry.types_state_b = entry.types_state_b, entry.types_state_a

    def drop_state_b(self, remove_dummies=False, atomname=None, resname=None, resid=None, atomtype=None):
        """
        Makes the topology non-alchemical again, just dropping
        all parameters for state B
        :param remove_dummies: bool, whether to remove A-state dummies
        :param atomname: str, name of the selected atom(s) for which state B will be dropped
        :param resname: str, name of the selected residue(s) for which state B will be dropped
        :param resid: int, number of the selected residue(s) for which state B will be dropped
        :param atomtype: str, type of the selected atom(s) for which state B will be dropped
        :return: None
        """
        if not remove_dummies:
            print("Warning: dropping all state B parameters, but keeping dummies (if exist). To remove all atoms with "
                  "names starting with D, rerun this fn with 'remove_dummies=True'.")
        if atomname or resname or resid or atomtype:
            selected = set()
            sub = self.get_subsection('atoms')
            for entries in [entry for entry in sub if isinstance(entry, gml.EntryAtom)]:
                criteria = all([(atomname is None or entries.atomname == atomname),
                                (resname is None or entries.resname == resname),
                                (resid is None or int(entries.resid) == int(resid)),
                                (atomtype is None or entries.type == atomtype)])
                if criteria:
                    selected.add(entries.num)
        else:
            selected = list(range(1, self.natoms + 1))
        for sub in self.subsections:
            for entry in sub:
                if (isinstance(entry, gml.EntryAtom) and entry.num in selected) \
                        or (isinstance(entry, gml.EntryBonded) and any([x in selected for x in entry.atom_numbers])):
                    if isinstance(entry, gml.EntryAtom) and entry.type_b is not None:
                        entry.type_b, entry.mass_b, entry.charge_b = 3 * [None]
                    elif isinstance(entry, gml.EntryBonded) and entry.params_state_b:
                        entry.params_state_b = []
                    if isinstance(entry, gml.EntryBonded) and entry.types_state_b is not None:
                        entry.types_state_b = None
        if remove_dummies:
            sub = self.get_subsection('atoms')
            dummies = [entry for entry in sub if isinstance(entry, gml.EntryAtom) and entry.atomname[0] == "D"]
            while dummies:
                to_remove = dummies[0]
                self.del_atom(to_remove.num)
                dummies = [entry for entry in sub if isinstance(entry, gml.EntryAtom) and entry.atomname[0] == "D"]
    
    def add_atom(self, atom_number, atom_name, atom_type, charge=0.0, resid=None, resname=None, mass=None, prnt=True):
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
            try:
                param_entry = [e for e in param_sect.get_subsection('atomtypes').entries
                               if isinstance(e, gml.EntryParam) and e.content[0] == atom_type][0]
                mass = param_entry.content[2]
            except IndexError:
                print("Could not assign mass for type {}, proceeding with 1.008 AU".format(atom_type))
                mass = 1.008
        fstring = subs_atoms.fstring
        if prnt:
            print(fstring.format(atom_number, atom_type, resid, resname, atom_name, atom_number, charge, mass).strip())
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
        """
        Removes an atom from the topology, as specified using
        topology numbering (1-based)
        :param atom_number: int, atom number in topology
        :param del_in_pdb: bool, whether to also remove in the bound PDB file
        :return: None
        """
        self._del_atom(atom_number)
        self._del_params(atom_number)
        self.offset_numbering(-1, atom_number)
        self.top.recalc_sys_params()
        if del_in_pdb:
            if self.top.pdb:
                for to_remove in self._match_pdb_to_top(atom_number):
                    self.top.pdb.delete_atom(to_remove)

    def swap_atom(self, atom_number, new_position, swap_in_pdb=True):
        """
        Changes the position of a chosen atom (1-based index atom_number)
        so that it now has index new_position (and other atoms are renumbered).
        If the topology has a corresponding structure, atoms can also be
        moved in the .pdb object.
        :param atom_number: int, atom to be moved (1-based)
        :param new_position: int, target index of the atom (1-based)
        :param swap_in_pdb: bool, whether to try moving the atom in Top.pdb
        :return: None
        """
        if swap_in_pdb:
            if self.top.pdb:
                if len(self._match_pdb_to_top(atom_number)) > 1:
                    raise RuntimeError("Two or more atoms in PDB matching the requested atom {} "
                                       "in .top".format(atom_number))
                elif len(self._match_pdb_to_top(atom_number)) == 0:
                    raise RuntimeError("Could not match .top atom {} to a corresponding PDB atom".format(atom_number))
                if len(self._match_pdb_to_top(new_position)) > 1:
                    raise RuntimeError("Two or more atoms in PDB matching the requested atom {} "
                                       "in .top".format(new_position))
                elif len(self._match_pdb_to_top(new_position)) == 0:
                    raise RuntimeError("Could not match .top atom {} to a corresponding PDB atom".format(new_position))
                old_loc = self._match_pdb_to_top(atom_number)[0]
                new_loc = self._match_pdb_to_top(new_position)
                atom = self.top.pdb.atoms.pop(old_loc-1)  # TODO check
                self.top.pdb.atoms.insert(new_loc + 1, atom)
        subsect_atoms = self.get_subsection('atoms')
        atom_entry_list = [e for e in subsect_atoms.entries]
        entry_ind = [n for n, e in enumerate(atom_entry_list) if isinstance(e, gml.EntryAtom)
                     and e.num == atom_number][0]
        self._hide_atom(atom_number, new_position)
        self.offset_numbering(-1, atom_number)
        self.offset_numbering(1, new_position)
        self._return_atom(new_position)
        entry_final_ind = [n for n, e in enumerate(atom_entry_list) if isinstance(e, gml.EntryAtom)][new_position - 1]
        entry = subsect_atoms.entries.pop(entry_ind)
        subsect_atoms.entries.insert(entry_final_ind, entry)

    def _hide_atom(self, old_pos, new_pos):
        subsect_atoms = self.get_subsection('atoms')
        chosen = [e for e in subsect_atoms.entries if isinstance(e, gml.EntryAtom) and e.num == old_pos][0]
        chosen.num = -new_pos
        for subs in ['bonds', 'angles', 'pairs', 'dihedrals', 'impropers', 'cmap']:
            try:
                subsection = self.get_subsection(subs)
                for entry in subsection:
                    if isinstance(entry, gml.EntryBonded):
                        if old_pos in entry.atom_numbers:
                            index = entry.atom_numbers.index(old_pos)
                            temp = list(entry.atom_numbers)
                            temp[index] = -new_pos
                            entry.atom_numbers = tuple(temp)
            except KeyError:
                pass

    def _return_atom(self, new_pos):
        subsect_atoms = self.get_subsection('atoms')
        chosen = [e for e in subsect_atoms.entries if isinstance(e, gml.EntryAtom) and e.num < 0][0]
        assert chosen.num == -new_pos
        chosen.num *= -1
        for subs in ['bonds', 'angles', 'pairs', 'dihedrals', 'impropers', 'cmap']:
            try:
                subsection = self.get_subsection(subs)
                for entry in subsection:
                    if isinstance(entry, gml.EntryBonded):
                        if any([x < 0 for x in entry.atom_numbers]):
                            if -new_pos in entry.atom_numbers:
                                index = entry.atom_numbers.index(-new_pos)
                                temp = list(entry.atom_numbers)
                                temp[index] *= -1
                                entry.atom_numbers = tuple(temp)
                            else:
                                print("Caution, found strange negative atom index in line {}".format(entry))
            except KeyError:
                pass
    
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
            self.top.recalc_sys_params()

    def merge_molecules(self, other):
        other.offset_numbering(self.natoms)
        self._merge_fields(other)
        self.top.sections.remove(other)
        # the stuff below works but is terribly ugly, we need to have API for manipulating content of Top.system
        system_setup = self.top.sections[-1].get_subsection('molecules')
        system_setup.entries = [e for e in system_setup if other.mol_name not in e]
        self.top.recalc_sys_params()

    def _merge_fields(self, other):
        # TODO important: watch for POSRES
        print('WARNING watch out for #ifdef POSRES keywords that might get misplaced')
        for subs in ['atoms', 'bonds', 'angles', 'pairs', 'dihedrals', 'impropers', 'cmap', 'position_restraints']:
            # TODO merge all subsections
            try:
                subsection_other = other.get_subsection(subs)
                subsection_own = self.get_subsection(subs)
                subsection_own.add_entries([entry for entry in subsection_other if entry])
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
        """
        Looks for FF parameters to be put for every bonded term in the topology,
        then adds them so that they can be explicitly seen/modified
        :param add_section: str, to which section should the FF params be added
        :return: None
        """
        if add_section == 'all':  # TODO optionally add type/atomname labels in comment
            subsections_to_add = ['bonds', 'angles', 'dihedrals', 'impropers']
        else:
            subsections_to_add = [add_section]
        for sub in subsections_to_add:
            try:
                subsections = [s for s in self.subsections if s.header == sub]
            except IndexError:
                pass
            else:
                for ssub in subsections:
                    ssub.add_ff_params()

    def find_used_ff_params(self, section='all'):
        used_params = []
        if section == 'all':  # TODO optionally add type/atomname labels in comment
            subsections_to_add = ['bonds', 'angles', 'dihedrals', 'impropers']
        else:
            subsections_to_add = [section]
        for sub in subsections_to_add:
            try:
                subsections = [s for s in self.subsections if s.header == sub]
            except IndexError:
                pass
            else:
                for ssub in subsections:
                    used_params.extend(ssub.find_used_ff_params())
        return used_params

    def find_missing_ff_params(self, add_section='all', fix_by_analogy=False, fix_B_from_A=False, fix_A_from_B=False):
        if add_section == 'all':  # TODO optionally add type/atomname labels in comment
            subsections_to_add = ['bonds', 'angles', 'dihedrals', 'impropers']
        else:
            subsections_to_add = [add_section]
        for sub in subsections_to_add:
            try:
                subsections = [s for s in self.subsections if s.header == sub]
            except IndexError:
                pass
            else:
                for ssub in subsections:
                    print(f"Searching in molecule {self.mol_name}, section {ssub}...")
                    ssub.find_missing_ff_params(fix_by_analogy, fix_B_from_A, fix_A_from_B)

    def label_types(self, add_section='all'):
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
                subsection.add_type_labels()

    def mutate_protein_residue(self, resid, target, rtp=None, mutate_in_pdb=True):
        alt_names = {('THR', 'OG'): 'OG1', ('THR', 'HG'): 'HG1', ('LEU', 'CD'): 'CD1', ('LEU', 'HD1'): 'HD11',
                     ('LEU', 'HD2'): 'HD12', ('LEU', 'HD3'): 'HD13', ('VAL', 'CG'): 'CG1', ('VAL', 'HG1'): 'HG11',
                     ('VAL', 'HG2'): 'HG12', ('VAL', 'HG3'): 'HG13',}
        if rtp is None and not self.top.rtp:
            print("Found the following .rtp files:\n")
            for n, i in enumerate(glob(self.top.gromacs_dir + '/*ff/[am][em]*rtp')):
                print('[', n + 1, '] ', i)
            rtpnum = input('\nPlease select one that contains the deserved charges and types:\n')
            try:
                rtpnum = int(rtpnum)
            except ValueError:
                raise RuntimeError('Not an integer: {}'.format(rtpnum))
            else:
                rtp = glob(self.top.gromacs_dir + '/*ff/[am][em]*rtp')[rtpnum-1]
        elif self.top.rtp:
            rtp = ''
        orig = self.atoms[self.select_atom('resid {} and name CA'.format(resid))]
        mutant = gml.ProteinMutant(orig.resname, target)
        targ = mutant.target_3l
        print("\n  Mutating residue {} (resid {}) into {}\n".format(orig.resname, resid, targ))
        atoms_add, hooks, _, _, extra_bonds, afters = mutant.atoms_to_add()
        atoms_remove = mutant.atoms_to_remove()
        types, charges, impropers, improper_type = self.parse_rtp(rtp)
        # some residue-specific modifications here
        if targ == 'HIS':
            targ = 'HSD' if ('HSD', 'CA') in types.keys() else 'HID'
        elif targ == 'GLY':
            self.atoms[self.select_atom('resid {} and name HA'.format(resid))].atomname = 'HA1'
        if orig.resname == 'GLY':
            self.atoms[self.select_atom('resid {} and name HA1'.format(resid))].atomname = 'HA'
        impropers_to_add = []
        impr_sub = self.get_subsection('impropers')
        atoms_sub = self.get_subsection('atoms')
        # first remove all unwanted atoms
        for at in atoms_remove:
            equivalents = {'OG': 'OG1', 'HG': 'HG1', 'HG1': 'HG11', 'HG2': 'HG12', 'HG3': 'HG13', 'CG': 'CG1',
                           'CD': 'CD1', 'HD': 'HD1', 'HD1': 'HD11', 'HD2': 'HD12', 'HD3': 'HD13'}
            print("Removing atom {} from resid {} in topology".format(at, resid))
            try:
                atnum = self.select_atom('resid {} and name {}'.format(resid, at))
            except RuntimeError:
                atnum = self.select_atom('resid {} and name {}'.format(resid, equivalents[at]))
            self.del_atom(self.atoms[atnum].num, del_in_pdb=False)
        for atom_add, hook, aft in zip(atoms_add, hooks, afters):
            print("Adding atom {} to resid {} in topology".format(atom_add, resid))
            # if there are ambiguities in naming (two or more options):
            if (targ, atom_add) in alt_names.keys():
                atom_add = alt_names[(targ, atom_add)]
            if isinstance(hook, tuple):
                for hk in hook:
                    try:
                        _ = self.select_atom('resid {} and name {}'.format(resid, hk))
                    except RuntimeError:
                        continue
                    else:
                        hook = hk
                        break
                else:
                    raise RuntimeError("Couldn't find any of the following atoms: {}".format(hook))
            hooksel = 'resid {} and name {}'.format(resid, hook)
            if isinstance(aft, tuple):
                for n, af in enumerate(aft):
                    try:
                        _ = self.select_atom('resid {} and name {}'.format(resid, af))
                    except RuntimeError:
                        continue
                    else:
                        aftnr = self.select_atom('resid {} and name {}'.format(resid, af))
                        break
                else:
                    raise RuntimeError("Couldn't find any of the following atoms: {}".format(aft))
            else:
                aftnr = self.select_atom('resid {} and name {}'.format(resid, aft))
            hnum = self.atoms[self.select_atom(hooksel)].num
            atnum = aftnr + 2
            # actual addition of atoms
            self.add_atom(atnum, atom_add, atom_type=types[(targ, atom_add)], charge=charges[(targ, atom_add)],
                          resid=orig.resid, resname=targ, mass=None, prnt=False)
            self.add_bond(hnum, atnum)
            for i in impropers[targ]:
                if atom_add in i and i not in impropers_to_add:
                    impropers_to_add.append(i)
        # changing resnames, charges and types according to .rtp
        for atom in self.select_atoms('resid {}'.format(resid)):
            self.atoms[atom].resname = targ
            try:
                self.atoms[atom].charge = charges[(targ, self.atoms[atom].atomname)]
            except KeyError:
                print("Couldn't find atom {} in RTP entry for residue {} - check charges and types "
                      "manually".format(self.atoms[atom].atomname, targ))
            self.atoms[atom].type = types[(targ, self.atoms[atom].atomname)]
        # bonds that close rings
        for bond in extra_bonds:
            xsel = 'resid {} and name {}'.format(resid, bond[0])
            ysel = 'resid {} and name {}'.format(resid, bond[1])
            xnum = self.atoms[self.select_atom(xsel)].num
            ynum = self.atoms[self.select_atom(ysel)].num
            self.add_bond(xnum, ynum)
        atoms_sub.get_dicts(force_update=True)
        # looking for new impropers
        for imp in impropers_to_add:
            if set(imp).intersection(set(atoms_add)):
                numbers = [atoms_sub.name_to_num[(resid, at)] for at in imp]
                new_str = '{:5d} {:5d} {:5d} {:5d} {:>5s}\n'.format(*numbers, improper_type)
                impr_sub.add_entry(gml.EntryBonded(new_str, impr_sub),
                                   position=1 + [n for n, e in enumerate(impr_sub) if isinstance(e, gml.EntryBonded)][-1])
        # repeating the mutation in the structure
        if mutate_in_pdb and self.top.pdb:
            pdb_atoms = self._match_pdb_to_top(self.atoms[self.select_atom('resid {} and name CA'.format(resid))].num)
            pdb_chains = [self.top.pdb.atoms[at].chain for at in pdb_atoms]
            if len(pdb_atoms) == 1:
                chain = '' if pdb_chains[0] == ' ' else pdb_chains[0]
                self.top.pdb.mutate_protein_residue(resid, target, chain)
            elif len(pdb_atoms) > 1:
                if any([pdb_chains[0] == pdb_chains[i] for i in range(1, len(pdb_chains))]):
                    print()
                    response = input("The topology entry {} corresponds to multiple entries in the PDB; should we add "
                                     "chains to PDB and retry? (y/n)\n".format(self.mol_name))
                    if response.lower() == 'y':
                        self.top.pdb.add_chains(maxwarn=-1)
                        pdb_chains = [self.top.pdb.atoms[at].chain for at in pdb_atoms]
                    else:
                        print("Mutated in .top, but not in .pdb; try running separately with "
                              "Pdb.mutate_protein_residue(), where chains can be specified separately")
                        return
                for ch in pdb_chains:
                    self.top.pdb.mutate_protein_residue(resid, target, ch)
        elif mutate_in_pdb and not self.top.pdb:
            print("No .pdb file bound to the topology, use Top.add_pdb() to add one")

    def parse_rtp(self, rtp):
        if self.top.rtp:
            return self.top.rtp['typedict'], self.top.rtp['chargedict'], self.top.rtp['impropers'], \
                   self.top.rtp['bondedtypes']
        chargedict, typedict = {}, {}
        impropers = {}
        bondedtypes = 0
        rtp_cont = [l for l in open(rtp) if not l.strip().startswith(';')]
        resname = None
        reading_atoms = False
        reading_impropers = False
        reading_bondedtypes = False
        for line in rtp_cont:
            if line.strip().startswith('[') and line.strip().split()[1] not in ['bondedtypes', 'atoms', 'bonds', 'impropers']:
                resname = line.strip().split()[1]
            if line.strip().startswith('[') and line.strip().split()[1] == 'atoms':
                reading_atoms = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'atoms':
                reading_atoms = False
            if line.strip().startswith('[') and line.strip().split()[1] == 'impropers':
                reading_impropers = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'impropers':
                reading_impropers = False
            if line.strip().startswith('[') and line.strip().split()[1] == 'bondedtypes':
                reading_bondedtypes = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'bondedtypes':
                reading_bondedtypes = False
            if len(line.strip().split()) > 3 and resname is not None and reading_atoms:
                typedict[(resname, line.strip().split()[0])] = line.strip().split()[1]
                chargedict[(resname, line.strip().split()[0])] = float(line.strip().split()[2])
            if len(line.strip().split()) > 3 and resname is not None and reading_impropers:
                if resname not in impropers.keys():
                    impropers[resname] = []
                impropers[resname].append(line.strip().split())
            if len(line.strip().split()) > 7 and resname is None and reading_bondedtypes:
                bondedtypes = line.strip().split()[3]
        # substitute CHARMM's HN for AMBER's H
        for k in list(typedict.keys()):
            if 'HN' in k:
                typedict[(k[0], 'H')] = typedict[k]
                chargedict[(k[0], 'H')] = chargedict[k]
            if 'HG1' in k:
                typedict[(k[0], 'HG')] = typedict[k]
                chargedict[(k[0], 'HG')] = chargedict[k]
        self.top.rtp['typedict'] = typedict
        self.top.rtp['chargedict'] = chargedict
        self.top.rtp['impropers'] = impropers
        self.top.rtp['bondedtypes'] = bondedtypes
        return self.top.rtp['typedict'], self.top.rtp['chargedict'], self.top.rtp['impropers'], \
            self.top.rtp['bondedtypes']

    def update_dicts(self): # TODO add whenever adding/removing molecules
        self.get_subsection('atoms').get_dicts(force_update=True)

    def list_bonds(self, by_types=False, by_params=False, returning=False):
        return self._list_bonded('bonds', by_types, by_params, returning)

    def list_angles(self, by_types=False, by_params=False, returning=False):
        return self._list_bonded('angles', by_types, by_params, returning)

    def list_impropers(self, by_types=False, by_params=False, returning=False):
        return self._list_bonded('impropers', by_types, by_params, returning)

    def list_dihedrals(self, by_types=False, by_params=False, returning=False):
        return self._list_bonded('dihedrals', by_types, by_params, returning)

    def _list_bonded(self, term, by_types, by_params, returning):
        self.update_dicts()
        subsection = self.get_subsection(term)
        returnable = []
        formatstring = {'bonds': "{:>5s} {:>5s}", 'angles': "{:>5s} {:>5s} {:>5s}",
                        'dihedrals': '{:>5s} {:>5s} {:>5s} {:>5s}', 'impropers': '{:>5s} {:>5s} {:>5s} {:>5s}'}
        for entry in subsection:
            if isinstance(entry, gml.EntryBonded):
                entry.read_types()
                if not by_params:
                    extra = ''
                    params = []
                else:
                    extra = '{:>12.5f} ' * len(entry.params_state_a)
                    params = entry.params_state_a
                if not returning:
                    if not by_types:
                        print((formatstring[term] + extra).format(*entry.atom_names, *params))
                    else:
                        print((formatstring[term] + extra).format(*entry.types_state_a, *params))
                else:
                    if not by_types:
                        returnable.append(entry.atom_names)
                    else:
                        returnable.append(entry.types_state_a)
        return None if not returning else returnable

    def alch_h_to_ch3(self, resid, orig_name, basename, ctype=None, htype=None, ccharge=None, hcharge=0.09,
                      dummy_type='DH', add_in_pdb=True):
        if ctype is None or htype is None:
            print("Which atomtypes should be used for the methyl group:\n")
            print("[ 1 ] CT/HC (Amber methyl)")
            print("[ 2 ] CT3/HA3 (Charmm methyl)")
            print("[ X/Y ] Use type X for carbon, type Y for hydrogen")
            sel = input("\n Please provide your selection:\n")
            if sel == '1':
                ctype, htype = 'CT', 'HC'
            elif sel == '2':
                ctype, htype = 'CT3', 'HA3'
            elif '/' in sel:
                ctype, htype = sel.split('/')
            else:
                raise RuntimeError("{} is not a valid selection".format(sel))
        self.add_dummy_def(dummy_type)
        orig = self.atoms[self.select_atoms('resid {}'.format(resid))[0]]
        ccharge = ccharge if ccharge is not None else round(orig.charge - 0.27, 4)
        atoms_add, hooks = [basename.replace('C', 'H') + str(i) for i in range(3)], 3 * [orig_name]
        for n, atom_add_hook in enumerate(zip(atoms_add, hooks), 1):
            atom_add, hook = atom_add_hook
            print("Adding atom {} to resid {} in the topology".format(atom_add, resid))
            hooksel = 'resid {} and name {}'.format(resid, orig_name)
            hnum = self.atoms[self.select_atom(hooksel)].num
            atnum = hnum + n
            self.add_atom(atnum, atom_add, atom_type=dummy_type, charge=0, resid=resid, resname=orig.resname, mass=1.008)
            self.add_bond(hnum, atnum)
            self.gen_state_b(atomname=atom_add, resid=resid, new_type=htype, new_charge=hcharge, new_mass=1.008)
        self.gen_state_b(atomname=orig_name, resid=resid, new_type=ctype, new_charge=ccharge, new_mass=12.0)
        if add_in_pdb and self.top.pdb:
            if len(self.top.system) > 1 or self.top.system[list(self.top.system.keys())[0]] > 1:
                raise RuntimeError("Adding groups in PDB only supported for systems containing one molecule")
            bonds = self.list_bonds(returning=True)
            hook = [j for i in bonds for j in i if orig_name in i and orig_name != j][0]
            aligns = [j for i in bonds for j in i if hook in i and hook != j and orig_name != j]
            aftnr = self.select_atom('resid {} and name {}'.format(resid, hook))
            for n, aliat in enumerate(zip(aligns, atoms_add), 1):
                ali, at = aliat
                self.top.pdb.insert_atom(aftnr+n, self.top.pdb.atoms[aftnr],
                                         atomsel='resid {} and name {}'.format(resid, at),
                                         hooksel='resid {} and name {}'.format(resid, orig_name), bondlength=1.1,
                                         p1_sel='resid {} and name {}'.format(resid, ali),
                                         p2_sel='resid {} and name {}'.format(resid, hook), atomname=at)

    def add_dummy_def(self, dummy_type):
        params = self.top.parameters
        atomtypes = params.get_subsection('atomtypes')
        dummy_entries = [e for e in atomtypes if isinstance(e, gml.EntryParam) and e.types[0] == dummy_type]
        if not dummy_entries:
            atomtypes.add_entry(gml.EntryParam('   {}     0        1.008  0.0000  A  0.000000000000  0.0000  '
                                               '\n'.format(dummy_type), atomtypes))

    def make_stateB_dummy(self, resid, orig_name, dummy_type='DH'):
        self.add_dummy_def(dummy_type)
        self.gen_state_b(atomname=orig_name, resid=resid, new_type=dummy_type, new_charge=0, new_mass=1.008)


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

    def sort_dihedrals(self):
        """
        Sorts dihedrals to make sure wildcards are
        moved to the very end of the file
        :return:
        """
        for sub in self.subsections:
            if 'dihedral' in sub.header:
                sub.sort()

    def clone_type(self, atomtype, prefix):
        """
        Generates an exact type of a selected atomtype,
        preserving all interactions with other types
        :param atomtype: str, atomtype to be duplicated
        :param prefix: str, new name will be generated as prefix + original atomtype
        :return: None
        """
        for sub in self.subsections:
            to_add = []
            for ent in sub:
                if isinstance(ent, gml.EntryParam) and atomtype in ent.types and prefix + atomtype not in ent.types:
                    to_add.append(ent)
            for entry in to_add:
                newlines = self.gen_clones(entry, atomtype, prefix)
                sub.add_entries([gml.EntryParam(line, sub) for line in newlines])
        self.sort_dihedrals()
        self._remove_symm_dupl(prefix)

    def clean_unused(self, used_params, section='all'):
        matchings = {'bonds': 'bondtypes', 'angles': 'angletypes', 'dihedrals': 'dihedraltypes',
                     'impropers': 'dihedraltypes'}
        if section == 'all':
            subs = list(matchings.values())
        else:
            subs = [matchings[section]]
        for sub in subs:
            ssects = [sb for sb in self.subsections if sb.header == sub]
            for ssect in ssects:
                new_entries = []
                for entry in ssect.entries:
                    if not isinstance(entry, gml.EntryParam) or entry.identifier in used_params:
                        new_entries.append(entry)
                ssect.entries = new_entries

    def _remove_symm_dupl(self, prefix):
        for sub in self.subsections:
            if 'dihedral' in sub.header:
                sub._remove_symm(prefix)

    def get_opt_dih(self):
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and int(sub.prmtype) == 9][0]
        return ss.get_opt_dih()

    def get_opt_dih_indices(self):
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and int(sub.prmtype) == 9][0]
        return ss.get_opt_dih_indices()

    def set_opt_dih(self, values):
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and int(sub.prmtype) == 9][0]
        ss.set_opt_dih(values)

    def add_nbfix(self, type1, type2, mod_sigma=0.0, mod_epsilon=0.0, action_default='x'):
        atp = self.get_subsection('atomtypes')
        sigma1, eps1, sigma2, eps2 = [None] * 4
        for entry in atp:
            if isinstance(entry, gml.EntryParam):
                if entry.types[0] == type1:
                    sigma1, eps1 = entry.params
                if entry.types[0] == type2:
                    sigma2, eps2 = entry.params
        if sigma1 is None:
            raise KeyError('Type {} was not found in the atomtype definitions'.format(type1))
        if sigma2 is None:
            raise KeyError('Type {} was not found in the atomtype definitions'.format(type2))
        new_sigma = 0.5*(sigma1 + sigma2) + mod_sigma
        new_epsilon = (eps1*eps2)**0.5 + mod_epsilon
        try:
            nbsub = self.get_subsection('nonbond_params')
        except KeyError:
            self.subsections.append(self._yield_sub(['[ nonbond_params ]']))
            nbsub = self.get_subsection('nonbond_params')
        comment = ''
        for entry in nbsub:
            if isinstance(entry, gml.EntryParam):
                if (entry.types[0], entry.types[1]) in [(type1, type2), (type2, type1)]:
                    action = action_default
                    while action not in 'mrt':
                        action = input("An entry already exists, shall we replace it (r), modify (m) or terminate (t)?")
                    if action == 't':
                        return
                    elif action == 'm':
                        new_sigma = entry.params[0] + mod_sigma
                        new_epsilon = entry.params[1] + mod_epsilon
                        comment = entry.comment
                    nbsub.remove_entry(entry)
        entry_line = "{} {} 1 {} {} ; sigma chg by {}, eps chg by {} {}".format(type1, type2, new_sigma, new_epsilon,
                                                                                mod_sigma, mod_epsilon, comment)
        nbsub.add_entry(gml.Subsection.yield_entry(nbsub, entry_line))

    @staticmethod
    def gen_clones(entry, atomtype, prefix):
        lines = []
        nchanges = entry.types.count(atomtype)
        changes = []
        for i in range(nchanges):
            changes.extend(SectionParam.gen_combs(nchanges, i + 1))
        for mod in changes:
            lines.append(SectionParam.mod_types(entry, mod, prefix, atomtype))
        return lines

    @staticmethod
    def gen_combs(count, tuples):
        return list(combinations(range(count), tuples))

    @staticmethod
    def mod_types(entry, mods, prefix, atomtype):
        line = str(entry)
        for num in mods[::-1]:
            indices = [i for i in range(len(line) - len(atomtype) + 1)
                       if line[i:i + len(atomtype)] == atomtype and (i == 0 or line[i-1].isspace())
                       and (i+len(atomtype) == len(line) or line[i+len(atomtype)].isspace())]
            line = line[:indices[num]] + prefix + line[indices[num]:]
        return line
