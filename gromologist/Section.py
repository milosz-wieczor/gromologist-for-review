import os
from itertools import product, combinations
from typing import Optional, Union
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
    
    def __init__(self, content: list, top: "gml.Top"):
        self.name = 'System'
        self.top = top
        self.conditional = 0
        self.dih_processed = False
        self.subsections = []
        content_split = self._split_content(content)
        excess_if = 0
        for cont in range(len(content_split)):
            self.subsections.append(self._yield_sub(content_split[cont]))
            if cont == 0:
                self.subsections[-1].conditional = self.conditional
            else:
                excess_if += self.count_ifs(content_split[cont-1]) - self.count_ifs(content_split[cont-1], endif=True)
                self.subsections[-1].conditional = self.conditional + excess_if
            if self.subsections[-1].conditional > 0:
                if isinstance(self.subsections[-1], gml.SubsectionParam):
                    self.top.print(f"Subection {str(self.subsections[-1])} recognized as conditional, will not be merged")
                elif isinstance(self.subsections[-1], gml.SubsectionBonded) \
                        or isinstance(self.subsections[-1], gml.SubsectionAtom):
                    self.top.print(f"Subection {str(self.subsections[-1])} "
                                   f"in molecule {str(self.subsections[-1].section.subsections[0].entries[1]).split()[0]} "
                                   f"recognized as conditional")
                else:
                    self.top.print(f"Subection {str(self.subsections[-1])} recognized as conditional")
    
    def __repr__(self) -> str:
        return "{} section with {} subsections".format(self.name, len(self.subsections))
    
    @staticmethod
    def _split_content(content: list) -> list:
        """
        Splits a block of text (list of strings passed to the __init__,
        corresponding to the entire content of the given section)
        into a list of blocs, each starting with a [ section_header ]
        :param content: list of strings, content of section
        :return: list of lists of strings, contents of individual subsections
        """
        special_lines = [n for n, l in enumerate(content) if l.strip().startswith('[')] + [len(content)]
        return [content[beg:end] for beg, end in zip(special_lines[:-1], special_lines[1:])]
        
    def _yield_sub(self, content: list):
        """
        A wrapper that will select which kind of subsection
        should be instantiated (generic, bonded, or params)
        :param content: list of strings, content of the subsection
        :return: a Subsection instance (or a derived class)
        """
        until = content[0].index(']')
        header = content[0][:until].strip().strip('[]').strip()
        if header in {'bonds', 'pairs', 'angles', 'dihedrals', 'settles', 'exclusions', 'cmap', 'position_restraints',
                      'virtual_sites2', 'constraints'}:
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

    @staticmethod
    def count_ifs(content: list, endif: bool = False) -> int:
        if endif:
            return len([ln for ln in content if ln.strip().startswith("#endif")])
        else:
            return len([ln for ln in content if ln.strip().startswith("#ifdef") or ln.strip().startswith("#ifndef")])

    def get_subsection(self, section_name: str):
        """
        Returns the specified subsection; we always need to run merge()
        on SectionParam first to avoid duplicates
        :param section_name: str, name of the subsection to be returned
        :return: gml.Subsection
        """
        ssect = [s for s in self.subsections if s.header == section_name]
        if len(ssect) == 0:
            raise KeyError("Subsection {} not found, check your molecule topology!".format(section_name))
        elif len(ssect) > 1:
            raise RuntimeError("Error: subsection {} duplicated in {}".format(section_name, str(self)))
        return ssect[0]

    def get_subsections(self, section_name: str) -> list:
        """
        Returns the list of specified subsections
        :param section_name: str
        :return: None
        """
        ssect = [s for s in self.subsections if s.header == section_name]
        if len(ssect) == 0:
            raise KeyError("Subsection {} not found, check your molecule topology!".format(section_name))
        return ssect


class SectionMol(Section):
    """
    This class should wrap the subsections of a single molecule
    (i.e. one [ moleculetype ], one [ atoms ], one [ bonds ] etc.)
    """
    
    def __init__(self, content_list: list, top: "gml.Top"):
        super().__init__(content_list, top)
        self.mol_name = self.get_subsection('moleculetype').molname
        self.name = '{} molecule'.format(self.mol_name)
        self._merge()
        
    def __repr__(self) -> str:
        return self.name

    @property
    def atoms(self) -> list:
        return self.get_subsection('atoms').entries_atom

    @property
    def natoms(self) -> int:
        return len(self.atoms)

    @property
    def nmols(self) -> int:
        return sum([mol_count[1] for mol_count in self.top.system if mol_count[0] == self.mol_name])

    @property
    def charge(self) -> float:
        return sum([a.charge for a in self.atoms])

    @property
    def mass(self) -> float:
        return sum([a.mass for a in self.atoms])

    @property
    def residues(self) -> list:
        resid = None
        reslist = []
        for at in self.atoms:
            if at.resid != resid:
                resid = at.resid
                reslist.append(f'{at.resname}-{at.resid}')
        return reslist

    def _merge(self):
        parsed_headers = set()
        to_remove = []
        for n, sub in enumerate(self.subsections):
            if sub.header not in parsed_headers:
                parsed_headers.add(sub.header)
            else:
                first_ssect = [s for s in self.subsections if s.header == sub.header][0]
                if first_ssect.conditional or sub.conditional:
                    print(f"Cannot merge sections {sub.header} because at least one is conditional, be careful")
                    continue
                first_ssect.add_entries([e for e in sub.entries if not e.is_header()])
                to_remove.append(n)
                first_ssect.prmtypes = first_ssect._check_parm_type()
        for ssect_index in to_remove[::-1]:
            self.subsections.pop(ssect_index)

    def _remove_subsection(self, header):
        """
        Removes a subsection with name 'header' (useful e.g. for removing position_restraints subsections)
        :param header: str, name of the subsection
        :return: None
        """
        # TODO if subs is conditional, remove the ifdefs???
        try:
            sub_ind = self.subsections.index(self.get_subsection(header))
        except KeyError:
            print(f'Could not find subsection {header}')
            return
        _ = self.subsections.pop(sub_ind)

    def set_type(self, type_to_set: str, atomname: str, resname: Optional[str] = None, resid: Optional[int] = None,
                 prefix: Optional[str] = None):
        """
        Sets a defined atomic type for all atoms fulfilling a criterion
        (residue name/id, atom name)
        :param type_to_set: str, the type to be set
        :param atomname: str, atoms with this atomname will have their type changed
        :param resname: str, atoms that will have their type changed will be restricted to this residue name
        :param resid: int, atoms that will have their type changed will be restricted to this residue id
        :param prefix: str, if specified then prefix + original type will be set instead of type_to_set
        :return: None
        """
        if resname is None:
            resnames = list({a.resname for a in self.atoms})
        else:
            if isinstance(resname, list) or isinstance(resname, tuple):
                resnames = resname
            else:
                resnames = [resname]
        if resid is None:
            resids = list({a.resid for a in self.atoms})
        else:
            if isinstance(resid, list) or isinstance(resid, tuple):
                resids = resid
            else:
                resids = [resid]
        for a in self.atoms:
            if a.resname in resnames and a.atomname == atomname and a.resid in resids:
                if prefix is None:
                    a.type = type_to_set
                else:
                    a.type = prefix + a.type

    def select_atoms(self, selection_string: str) -> list:
        """
        Returns atoms' indices according to the specified selection string
        :param selection_string: str, a VMD-compatible selection
        :return: list, 0-based indices of atoms compatible with the selection
        """
        sel = gml.SelectionParser(self)
        return sel(selection_string)

    def select_atom(self, selection_string: str) -> int:
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

    def get_atoms(self, selection_string: str) -> list:
        return [self.atoms[i] for i in self.select_atoms(selection_string)]

    def get_atom(self, selection_string: str) -> "gml.EntryAtom":
        return self.atoms[self.select_atom(selection_string)]

    def print_molecule(self):
        sub = self.get_subsection('atoms')
        for entry in sub:
            print(str(entry), end='')

    @property
    def is_alchemical(self) -> bool:
        sect = self.get_subsection('atoms')
        for ent in sect.entries_atom:
            if ent.type_b is not None:
                return True
        return False

    @property
    def is_water(self) -> bool:
        return bool(self.get_atoms('water'))

    @property
    def is_protein(self) -> bool:
        return bool(self.get_atoms('protein'))

    @property
    def is_nucleic(self) -> bool:
        return bool(self.get_atoms('nucleic'))

    def _patch_alch(self):
        self.update_dicts()
        resnames = {a.resname for a in self.atoms}
        if 'DTS' in resnames:
            selected = [x+1 for x in self.select_atoms('resname DTS')]
            sect = self.get_subsection('bonds')
            for entry in sect.entries_bonded:
                if all([x in selected for x in entry.atom_numbers]):
                    entry.read_types()
                    if entry.atom_names == ("C5'", "O5'") or entry.atom_names == ("O5'", "C5'"):
                        if 'OS' in entry.types_state_a:
                            entry.params_state_b = [0.40000, 0.0]
                            entry.comment = entry.comment.rstrip() + ' fixed\n' if entry.comment else '; fixed\n'
                        else:
                            raise RuntimeError("Expected type OS in state A for bond C5'-O5' but found types {}, "
                                               "aborting".format(' '.join(entry.types_state_a)))
            sect = self.get_subsection('pairs')
            to_remove = []
            for n, entry in enumerate(sect):
                if isinstance(entry, gml.EntryBonded) and all([x in selected for x in entry.atom_numbers]):
                    entry.read_types()
                    if entry.atom_names == ("P", "Ox5'") or entry.atom_names == ("P", "DO5'") or \
                            entry.atom_names == ("O5'", "Hx5'") or entry.atom_names == ("O5'", "DH5'"):
                        to_remove.append(n)
            for en in to_remove[::-1]:
                _ = sect.entries.pop(en)
            atoms = sorted([x+1 for x in self.select_atoms("resname DTS and name O5' C5'")])
            assert len(atoms) % 2 == 0
            for pair in range(len(atoms)//2):
                self.nullify_bonded(*atoms[2*pair:2*pair+2], 'angles')
                self.nullify_bonded(*atoms[2 * pair:2 * pair + 2], 'dihedrals')
        if 'DTD' in resnames or 'DTE' in resnames:
            selected = [x+1 for x in self.select_atoms('resname DTD DTE')]
            sect = self.get_subsection('bonds')
            for entry in sect.entries_bonded:
                if all([x in selected for x in entry.atom_numbers]):
                    entry.read_types()
                    if entry.types_state_a == ("DUM_CT", "DUM_C2") \
                            and entry.atom_numbers[1] - entry.atom_numbers[0] > 30:
                        entry.params_state_a[1] = 0.0
                        entry.comment = entry.comment.rstrip() + ' fixed\n' if entry.comment else '; fixed\n'
        if 'DTX' in resnames or 'DTY' in resnames:
            selected = [x+1 for x in self.select_atoms('resname DTX DTY')]
            sect = self.get_subsection('bonds')
            for entry in sect.entries_bonded:
                if all([x in selected for x in entry.atom_numbers]):
                    entry.read_types()
                    if entry.types_state_a == ("DUM_CT", "DUM_CT") \
                            and entry.atom_numbers[1] - entry.atom_numbers[0] > 30:
                        entry.params_state_a[1] = 0.0
                        entry.comment = entry.comment.rstrip() + ' fixed\n' if entry.comment else '; fixed\n'

    def nullify_bonded(self, atom1, atom2, subsection, stateB=True):
        subs = self.get_subsection(subsection)
        for entry in subs.entries_bonded:
            if atom1 in entry.atom_numbers and atom2 in entry.atom_numbers:
                entry.read_types()
                if subsection == 'angles' and any(x.startswith('D') for x in entry.types_state_a):
                    continue
                if not entry.params_state_b:
                    entry.params_state_b = entry.params_state_a[:]
                if stateB:
                    entry.params_state_b[1] = 0.0
                else:
                    entry.params_state_a[1] = 0.0
                entry.comment = entry.comment.rstrip() + ' fixed\n' if entry.comment else '; fixed\n'

    def offset_numbering(self, offset: int, startfrom: int = 0):
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

    def _offset_atoms(self, offset: int, startfrom: int):
        """
        Offsets atoms in the [ atoms ] section
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        subsection = self.get_subsection('atoms')
        for entry in subsection.entries_atom:
            if entry.num >= startfrom:
                entry.num += offset

    def _offset_params(self, offset: int, startfrom: int):
        """
        Offsets atomic numbering in all parameter sections,
        e.g., [ bonds ]
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        for sub_name in [s.header for s in self.subsections if s.header != 'atoms']:
            subsection = self.get_subsection(sub_name)
            try:
                for entry in subsection.entries_bonded:
                    entry.atom_numbers = tuple(n + (offset * (n >= startfrom)) for n in entry.atom_numbers)
            except AttributeError:
                continue

    def add_disulfide(self, resid1: int, resid2: int, other: Optional["gml.SectionMol"] = None, rtp: Optional[str] = None):
        """
        Adds a bond between the SG atoms of a selected residue pair (either within
        or between molecules), removes the extra hydrogens, and adjusts atom types
        or charges according to the .rtp file chosen
        :param resid1: int, the resid of 1st residue to be linked
        :param resid2: int, the resid of 2st residue to be linked (can be from a different molecule)
        :param other: SectionMol (molecule), optional; if specified, the disulfide will be intermolecular
        :param rtp: str, optional; if specified, the given .rtp will be used (otherwise interactive)
        :return: None
        """
        other = self if other is None else other
        s1 = self.get_atom(f'resid {resid1} and name SG')
        h1 = self.get_atom(f'resid {resid1} and name HG HG1')
        s2 = other.get_atom(f'resid {resid2} and name SG')
        h2 = other.get_atom(f'resid {resid2} and name HG HG1')
        assert s1.resname == s2.resname == 'CYS'
        self.del_atom(h1.num)
        other.del_atom(h2.num)
        found = self.find_rtp(rtp)
        types, charges, dihedrals, impropers, improper_type = self.parse_rtp(found)
        disulf_names = {'CYX', 'CYS2'}
        disulf = set([x[0] for x in types.keys()]).intersection(disulf_names)
        if len(disulf) == 0:
            self.top.rtp = {}
            raise RuntimeError(f"None of the residues {disulf_names} found in the .rtp file {found}")
        disulf_resname = disulf.pop()
        for mol, res in zip([self, other], [resid1, resid2]):
            for atom in mol.get_atoms(f'resid {res}'):
                atom.type = types[(disulf_resname, atom.atomname)]
                atom.charge = charges[(disulf_resname, atom.atomname)]
                atom.resname = disulf_resname
        self.merge_two(other, s1.num, s2.num)

    def add_coordinated_ion(self, resid1: int, resid2: int, other: Optional["gml.SectionMol"] = None, rtp: Optional[str] = None):
        """
        Adds a bond between the cysteine sulfur or histidine nitrogen and (usually)
        a transition metal ion like Zn or Fe, and adjusts atom types
        or charges according to the .rtp file chosen
        :param resid1: int, the resid of 1st residue to be linked (amino acid)
        :param resid2: int, the resid of 2st residue to be linked (ion/heme; can be from a different molecule)
        :param other: SectionMol (molecule), optional; if specified, the bond will be intermolecular
        :param rtp: str, optional; if specified, the given .rtp will be used (otherwise interactive)
        :return: None
        """
        other = self if other is None else other
        c1 = self.get_atom(f'resid {resid1} and name CA')
        if c1.resname not in {'CYS', 'HIS', 'HID', 'HIE', 'HSD', 'HSE'}:
            raise RuntimeError(f"Residue 1 has to be cysteine or histidine, not {c1.resname}")
        x2 = other.get_atom(f'resid {resid2} and name ZN ZN2 FE')
        if x2.resname not in {'ZN', 'ZN2', 'FE', 'HEM', 'HEME'}:
            raise RuntimeError(f"Residue 2 has to be a zinc or iron/heme, not {x2.resname}")
        if c1.resname == 'CYS':
            found = self.find_rtp(rtp)
            types, charges, dihedrals, impropers, improper_type = self.parse_rtp(found)
            h1 = self.get_atom(f'resid {resid1} and name HG HG1')
            x1 = self.get_atom(f'resid {resid1} and name SG')
            self.del_atom(h1.num)
            sulf_names = {'CYM'}
            sulf = set([x[0] for x in types.keys()]).intersection(sulf_names)
            if len(sulf) == 0:
                self.top.rtp = {}
                raise RuntimeError(f"None of the residues {sulf_names} found in the .rtp file {found}")
            sulf_resname = sulf.pop()
            for atom in self.get_atoms(f'resid {resid1}'):
                atom.type = types[(sulf_resname, atom.atomname)]
                atom.charge = charges[(sulf_resname, atom.atomname)]
                atom.resname = sulf_resname
        else:
            try:
                self.get_atom(f'resid {resid1} and name HE2')
            except RuntimeError:
                try:
                    self.get_atom(f'resid {resid1} and name HD1')
                except RuntimeError:
                    raise RuntimeError(f"In residue {c1.resname} {resid1}, neither HE2 nor HD1 were found, can't "
                                       f"identify histidine protonation state")
                else:
                    x1 = self.get_atom(f'resid {resid1} and name NE2')
            else:
                x1 = self.get_atom(f'resid {resid1} and name ND1')
        self.merge_two(other, x1.num, x2.num)

    def gen_state_b(self, atomname: Optional[str] = None, resname: Optional[str] = None,
                    resid: Optional[int] = None, atomtype: Optional[str] = None, new_type: Optional[str] = None,
                    new_charge: Optional[float] = None, new_mass: Optional[float] = None):
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
        for entry in sub.entries_atom:
            criteria = all([(atomname is None or entry.atomname == atomname),
                            (resname is None or entry.resname == resname),
                            (resid is None or int(entry.resid) == int(resid)),
                            (atomtype is None or entry.type == atomtype)])
            if criteria:
                entry.type_b = new_type if new_type is not None else entry.type
                entry.mass_b = new_mass if new_mass is not None else entry.mass
                entry.charge_b = new_charge if new_charge is not None else entry.charge
        self.update_dicts()

    def drop_state_a(self, remove_dummies: bool = False, atomname: Optional[str] = None, resname: Optional[str] = None,
                     resid: Optional[int] = None, atomtype: Optional[str] = None):
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
            for entry in sub.entries_atom:
                criteria = all([(atomname is None or entry.atomname == atomname),
                                (resname is None or entry.resname == resname),
                                (resid is None or int(entry.resid) == int(resid)),
                                (atomtype is None or entry.type == atomtype)])
                if criteria:
                    selected.add(entry.num)
        else:
            selected = list(range(1, self.natoms+1))
        if remove_dummies:
            sub = self.get_subsection('atoms')
            dummies = [entry for entry in sub.entries_atom if entry.type_b and entry.type_b[0] == "D"
                       and entry.num in selected]
            while dummies:
                to_remove = dummies[-1]
                self.del_atom(to_remove.num)
                dummies = [entry for entry in sub.entries_atom if entry.type_b and entry.type_b[0] == "D"
                           and entry.num in selected]
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
        self.update_dicts()

    def swap_states(self, atomname: Optional[str] = None, resname: Optional[str] = None,
                     resid: Optional[int] = None, atomtype: Optional[str] = None):
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
            for entry in sub.entries_atom:
                criteria = all([(atomname is None or entry.atomname == atomname),
                                (resname is None or entry.resname == resname),
                                (resid is None or int(entry.resid) == int(resid)),
                                (atomtype is None or entry.type == atomtype)])
                if criteria:
                    selected.add(entry.num)
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
        self.update_dicts()

    def drop_state_b(self, remove_dummies: bool = False, atomname: Optional[str] = None, resname: Optional[str] = None,
                     resid: Optional[int] = None, atomtype: Optional[str] = None):
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
            for entry in sub.entries_atom:
                criteria = all([(atomname is None or entry.atomname == atomname),
                                (resname is None or entry.resname == resname),
                                (resid is None or int(entry.resid) == int(resid)),
                                (atomtype is None or entry.type == atomtype)])
                if criteria:
                    selected.add(entry.num)
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
            dummies = [entry for entry in sub.entries_atom if entry.type[0] == "D" and entry.num in selected]
            while dummies:
                to_remove = dummies[-1]
                self.del_atom(to_remove.num)
                dummies = [entry for entry in sub.entries_atom if entry.type[0] == "D" and entry.num in selected]
        self.update_dicts()
    
    def add_atom(self, atom_number: int, atom_name: str, atom_type: str, charge: float = 0.0,
                 resid: Optional[int] = None, resname: Optional[str] = None, mass: Optional[float] = None,
                 prnt: bool = True):
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
        :param prnt: bool, whether to print the atom being added
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
        if mass is None:
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
        self.update_dicts()
    
    def del_atom(self, atom_number: int, del_in_pdb: bool = True, renumber_in_pdb: bool = True):
        """
        Removes an atom from the topology, as specified using
        topology numbering (1-based)
        :param atom_number: int, atom number in topology
        :param del_in_pdb: bool, whether to also remove in the bound PDB file
        :param renumber_in_pdb: bool, whether to renumber atoms in the bound PDB file
        :return: None
        """
        if atom_number > self.natoms:
            raise RuntimeError(f"Can't remove atom {atom_number}, molecule only has {self.natoms} atoms")
        if self.top.pdb:
            matched = self._match_pdb_to_top(atom_number) if del_in_pdb else []
        self._del_atom(atom_number)
        self._del_params(atom_number)
        self.offset_numbering(-1, atom_number)
        self.update_dicts()
        # checking correct numbers, can add more checks in the future
        self._check_correct()
        if del_in_pdb:
            if self.top.pdb:
                for to_remove in matched:
                    self.top.pdb.delete_atom(to_remove)
                if renumber_in_pdb:
                    self.top.pdb.renumber_atoms()

    def swap_atom(self, atom_number: int, new_position: int, swap_in_pdb: bool = True):
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
                self.top.pdb.renumber_atoms()
                old_locs = self._match_pdb_to_top(atom_number)
                new_locs = self._match_pdb_to_top(new_position)
                for old_loc, new_loc in zip(old_locs[::-1], new_locs[::-1]):
                    atom = self.top.pdb.atoms.pop(old_loc - 1)
                    self.top.pdb.atoms.insert(new_loc - 1, atom)
                self.top.pdb.renumber_atoms()
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
        self.update_dicts()

    def _hide_atom(self, old_pos: int, new_pos: int):
        subsect_atoms = self.get_subsection('atoms')
        chosen = [e for e in subsect_atoms.entries if isinstance(e, gml.EntryAtom) and e.num == old_pos][0]
        chosen.num = -new_pos
        for subs in ['bonds', 'angles', 'pairs', 'dihedrals', 'cmap']:
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

    def _return_atom(self, new_pos: int):
        subsect_atoms = self.get_subsection('atoms')
        chosen = [e for e in subsect_atoms.entries if isinstance(e, gml.EntryAtom) and e.num < 0][0]
        assert chosen.num == -new_pos
        chosen.num *= -1
        for subs in ['bonds', 'angles', 'pairs', 'dihedrals', 'cmap']:
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
    
    def _match_pdb_to_top(self, atom_number: int) -> list:
        """
        Returns a list of PDB atom indices (assuming .top matches .pdb)
        that correspond to the specified atom_number in the molecule topology
        :param atom_number: int, atom number in self (1-based)
        :return: list, PDB atom serials (1-based)
        """
        if not self.top.pdb:
            raise ValueError("No PDB object matched to the currently processed topology")
        count = atom_number - 1
        pdb_atom_serials = []
        for mol_count in self.top.system:
            for _ in range(mol_count[1]):
                if mol_count[0] == self.mol_name:
                    pdb_atom_serials.append(self.top.pdb.atoms[count].serial)
                count += self.top.get_molecule(mol_count[0]).natoms
        if len(pdb_atom_serials) > self.top.nmol(self.mol_name):
            raise RuntimeError("Too many atoms atoms in PDB matching the requested atom {} "
                               "in .top".format(atom_number))
        elif len(pdb_atom_serials) == 0:
            raise RuntimeError("Could not match .top atom {} to a corresponding PDB atom".format(atom_number))
        return pdb_atom_serials
        
    def _del_atom(self, atom_number: int):
        subsect_atoms = self.get_subsection('atoms')
        chosen = [e for e in subsect_atoms.entries if isinstance(e, gml.EntryAtom) and e.num == atom_number][0]
        subsect_atoms.entries.remove(chosen)
    
    def _del_params(self, atom_number: int):
        for subs in ['bonds', 'angles', 'pairs', 'dihedrals', 'cmap', 'position_restraints']:
            try:
                subsection = self.get_subsection(subs)
                to_del = []
                for entry in subsection.entries_bonded:
                    if atom_number in entry.atom_numbers:
                        to_del.append(entry)
                for entry in to_del:
                    subsection.entries.remove(entry)
            except KeyError:
                pass

    def _check_correct(self):
        for subs in ['bonds', 'angles', 'pairs', 'dihedrals', 'cmap', 'position_restraints']:
            try:
                subsection = self.get_subsection(subs)
                for entry in subsection.entries_bonded:
                    if any([e > self.natoms for e in entry.atom_numbers]):
                        raise RuntimeError(f"Entry {entry} is invalid, only {self.natoms} atoms in the system")
            except KeyError:
                pass

    @property
    def bonds(self):
        subsection = self.get_subsection('bonds')
        return [entry.atom_numbers for entry in subsection.entries_bonded]

    def add_bond(self, first_atom: int, second_atom: int):
        """
        This is just an alias for merge_two if bond is intramolecular
        """
        self.merge_two(self, first_atom, second_atom)

    def select_atoms_bonded_to(self, serial: int) -> list:
        """
        Returns indices of atoms that are bonded to atom with a specified number
        :param serial: int, 1-based index of the central atom
        :return: list of int, 1-based indices of the atoms bonded to it
        """
        bonded = []
        for bond in self.bonds:
            if serial in bond:
                bonded.append([x for x in bond if x != serial][0])
        return bonded

    def get_atoms_bonded_to(self, serial: int) -> list:
        """
        As self.select_atoms_bonded_to(), but instead of indices returns atom instances
        :param serial: int, 1-based index of the central atom
        :return: list of gml.EntryAtom, representations of the atoms bonded to it
        """
        return [self.atoms[i-1] for i in self.select_atoms_bonded_to(serial)]

    def merge_two(self, other: "gml.SectionMol", anchor_own: int, anchor_other: int):
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
        if self.top.pdb:
            print("WARNING: if merging molecules that are not consecutive, make sure to adjust your PDB numbering")

    def merge_molecules(self, other: "gml.SectionMol"):
        other.offset_numbering(self.natoms)
        self._merge_fields(other)
        self.top.sections.remove(other)
        # the stuff below works but is terribly ugly, we need to have API for manipulating content of Top.system
        system_setup = self.top.sections[-1].get_subsection('molecules')
        system_setup.entries = [e for e in system_setup if other.mol_name not in e]

    def _merge_fields(self, other: "gml.SectionMol"):
        self.top.print('WARNING watch out for #ifdef POSRES keywords that might get misplaced')
        for subs in ['atoms', 'bonds', 'angles', 'pairs', 'dihedrals', 'cmap', 'position_restraints']:
            # TODO merge all subsections
            # TODO check for a "conditional" attribute of a section
            try:
                subsection_other = other.get_subsection(subs)
                subsection_own = self.get_subsection(subs)
                subsection_own.add_entries([entry for entry in subsection_other if entry])
                subsection_other.section = subsection_own.section
                self.top.print(f"Merging sections {subs} from two molecules")
            except KeyError:
                pass
    
    def _make_bond(self, atom_own: int, atom_other: int, other: "gml.SectionMol"):
        try:
            _ = other.bonds
        except KeyError:
            other.subsections.append(self._yield_sub([f"[ bonds ]\n"]))
        new_bond = [tuple(sorted([int(atom_own), int(atom_other)]))]
        new_angles = self._generate_angles(other, atom_own, atom_other)
        new_pairs, new_dihedrals = self._generate_14(other, atom_own, atom_other)
        # TODO remove overlapping pairs between new_bond/new_angles and new_pairs for 4- and 5-membered rings
        # or do we really need it?
        for sub, entries in zip(['bonds', 'pairs', 'angles', 'dihedrals'],
                                [new_bond, new_pairs, new_angles, new_dihedrals]):
            subsection = self.get_subsection(sub)
            subsection.add_entries([gml.EntryBonded(subsection.fstring.format(*entry, subsection.prmtypes[0]), subsection)
                                    for entry in entries])
            try:
                other.get_subsection(sub)
            except KeyError:
                other.subsections.append(self._yield_sub([f"[ {sub} ]\n"]))

    def remove_bond(self, at1: int, at2: int):
        bond_to_remove = [(at1, at2)]
        if not (bond_to_remove[0] in self.bonds or tuple(x for x in bond_to_remove[0][::-1]) in self.bonds):
            raise RuntimeError("Bond between atoms {} and {} not found in the topology".format(at1, at2))
        angles_to_remove = self._generate_angles(self, at1, at2)
        pairs_to_remove, dihedrals_to_remove = self._generate_14(self, at1, at2)
        impropers = self.get_subsection('dihedrals')
        impropers_to_remove = []
        for n, entry in enumerate(impropers.entries):
            if isinstance(entry, gml.EntryBonded) and at1 in entry.atom_numbers and at2 in entry.atom_numbers:
                impropers_to_remove.append(n)
        for n in impropers_to_remove[::-1]:
            _ = impropers.entries.pop(n)
        try:
            cmaps = self.get_subsection('cmap')
        except KeyError:
            pass
        else:
            cmap_to_remove = []
            for n, entry in enumerate(cmaps.entries):
                if isinstance(entry, gml.EntryBonded) and at1 in entry.atom_numbers and at2 in entry.atom_numbers:
                    cmap_to_remove.append(n)
            for n in cmap_to_remove[::-1]:
                _ = cmaps.entries.pop(n)

        def match(seq1, seqlist):
            for seq2 in seqlist:
                if all(i == j for i, j in zip(seq1, seq2)) or all(i == j for i, j in zip(seq1, seq2[::-1])):
                    return True
            return False

        for sub, removable in zip(['bonds', 'pairs', 'angles', 'dihedrals'],
                                  [bond_to_remove, pairs_to_remove, angles_to_remove, dihedrals_to_remove]):
            subsection = self.get_subsection(sub)
            to_remove = []
            for n, e in enumerate(subsection.entries):
                if isinstance(e, gml.EntryBonded) and match(e.atom_numbers, removable):
                    to_remove.append(n)
            for n in to_remove[::-1]:
                _ = subsection.entries.pop(n)

    def _generate_angles(self, other: "gml.SectionMol", atom_own: int, atom_other: int):
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

    def _generate_14(self, other: "gml.SectionMol", atom_own: int, atom_other: int) -> (list, list):
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
        # cleanup same-1-4 cases for 3-membered rings etc
        # TODO also remove pairs for 4- and 5-membered rings
        # TODO make sure there are no replicates wrt existing terms
        new_pairs = [pair for pair in new_pairs if pair[0] != pair[-1]]
        new_dihedrals = [dih for dih in new_dihedrals if dih[0] != dih[-1]]
        return new_pairs, new_dihedrals
    
    def add_ff_params(self, add_section: str = 'all', force_all: bool = False):
        """
        Looks for FF parameters to be put for every bonded term in the topology,
        then adds them so that they can be explicitly seen/modified
        :param add_section: str, to which section should the FF params be added
        :param force_all: bool, whether to overwrite existing parameters
        :return: None
        """
        if add_section == 'all':
            subsections_to_add = ['bonds', 'angles', 'dihedrals']
        else:
            subsections_to_add = [add_section]
        for sub in subsections_to_add:
            try:
                subsections = [s for s in self.subsections if s.header == sub]
            except IndexError:
                pass
            else:
                for ssub in subsections:
                    ssub.add_ff_params(force_all)

    def find_used_ff_params(self, section: str = 'all') -> list:
        used_params = []
        if section == 'all':
            subsections_to_add = ['bonds', 'angles', 'dihedrals', 'cmap']
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

    def find_missing_ff_params(self, add_section: str = 'all', fix_by_analogy: bool = False, fix_B_from_A: bool = False,
                               fix_A_from_B: bool = False, fix_dummy: bool = False, once: bool = False):
        if add_section == 'all':
            subsections_to_add = ['bonds', 'angles', 'dihedrals']
        else:
            subsections_to_add = [add_section]
        self.printed = []
        self.get_subsection('atoms').check_defined_types()
        for sub in subsections_to_add:
            try:
                subsections = [s for s in self.subsections if s.header == sub]
            except IndexError:
                pass
            else:
                for ssub in subsections:
                    print(f"Searching in molecule {self.mol_name}, section {ssub}...")
                    ssub.find_missing_ff_params(fix_by_analogy, fix_B_from_A, fix_A_from_B, fix_dummy, once)
        del self.printed

    def label_types(self, add_section: str = 'all'):
        if add_section == 'all':
            subsections_to_add = ['bonds', 'angles', 'dihedrals']
        else:
            subsections_to_add = [add_section]
        for sub in subsections_to_add:
            try:
                subsection = [s for s in self.subsections if s.header == sub][0]
            except IndexError:
                pass
            else:
                subsection.add_type_labels()

    def hydrogen_mass_repartitioning(self, hmass: float = 4.032):
        """
        Repartitions the masses from heavy atoms to hydrogens to ensure that each
        hydrogen has the desired mass; this enables the use of a 4-fs time step in
        standard MD simulations
        :param hmass: float, desired mass of the hydrogen atom; default is 4.032 (as set by -heavyh in pdb2gmx)
        :return: None
        """
        for atom in self.atoms:
            if not atom.ish:
                hydrogens = [at for at in self.get_atoms_bonded_to(atom.num) if at.ish]
                hmasses_diff = [hmass - a.mass for a in hydrogens]
                atom.mass -= round(sum(hmasses_diff), 6)
                assert atom.mass > 2, f"After repartitioning atom {atom} has mass lower than 2 Daltons which beats " \
                                      f"the purpose of repartitioning"
                for hd, ha in zip(hmasses_diff, hydrogens):
                    ha.mass += hd

    def add_posres(self, keyword: Union[str, None] = 'POSRES', value: float = 500.0, selection=None):
        """
        Adds a position restraint section to the topology
        :param keyword: conditional keyword that will be used in the #ifdef directive, default is POSRES;
        if None, the POSRES entry won't be conditional
        :param value: value of the force constant, default is 500
        :param selection: position_restrains will be added for atoms that fit this selection (e.g. "name CA"); note
        that hydrogens are excluded anyway
        :return: None
        """
        try:
            _ = self.get_subsection('position_restraints')
        except KeyError:
            try:
                _ = value[0]
            except TypeError:
                value = 3 * [value]
            content = ['[ position_restraints ]', f'#ifdef {keyword}' if keyword is not None else '',
                       '; ai  funct  fcx    fcy    fcz']
            atoms = self.atoms if selection is None else self.get_atoms(selection)
            for atom in atoms:
                if len(atom.atomname) > 1:
                    if not atom.ish:
                        content.append(f"{atom.num:5}    1 {value[0]:5} {value[1]:5} {value[2]:5}")
            content.append("#endif\n" if keyword is not None else '')
            self.subsections.append(gml.SubsectionBonded(content, self))
        else:
            self.top.print(f"[ position_restraints ] already present in molecule {self.mol_name}, skipping")

    def find_rtp(self, rtp: Union[str, None]) -> str:
        """
        Looks for aminoacids.rtp or merged.rtp in local files (*ff/*rtp or *rtp)
        and in the Gromacs directory, then allows to interactively choose which one to use
        :param rtp: str, path to the .rtp file (if applicable)
        :return: str, path to the .rtp file found
        """
        if rtp is None and not self.top.rtp:
            found = glob(self.top.gromacs_dir + '/*ff/[am][em]*rtp') + glob(os.getcwd() + '/*ff/[am][em]*rtp') + \
                    glob(os.getcwd() + '/[am][em]*rtp')
            if not found:
                raise RuntimeError("No .rtp files found locally or in default Gromacs dirs. Please set "
                                   "rtp=/path/to/rtp/file")
            print("Found the following .rtp files:\n")
            for n, i in enumerate(found):
                print('[', n + 1, '] ', i)
            rtpnum = input('\nPlease select one that contains the deserved charges and types:\n')
            try:
                rtpnum = int(rtpnum)
            except ValueError:
                raise RuntimeError('Not an integer: {}'.format(rtpnum))
            else:
                rtp = found[rtpnum-1]
        elif self.top.rtp:
            rtp = None
        return rtp

    def mutate_protein_residue(self, resid: int, target: str, rtp: Optional[str] = None, mutate_in_pdb: bool = True):
        """
        Mutates an amino acid to a different one, optionally in the topology
        and structure simultaneously
        :param resid: int, number of the residue to be mutated
        :param target: str, a single-letter code of the new residue to be introduced
        :param rtp: str, path to the .rtp file that will be used to read atom properties
        :param mutate_in_pdb: bool, whether to attempt modifying the same residue in the associated Pdb
        :return: None
        """
        alt_names = {('THR', 'OG'): 'OG1', ('THR', 'HG'): 'HG1', ('LEU', 'CD'): 'CD1', ('LEU', 'HD1'): 'HD11',
                     ('LEU', 'HD2'): 'HD12', ('LEU', 'HD3'): 'HD13', ('VAL', 'CG'): 'CG1', ('VAL', 'HG1'): 'HG11',
                     ('VAL', 'HG2'): 'HG12', ('VAL', 'HG3'): 'HG13'}
        rtp = self.find_rtp(rtp) if rtp is None else rtp
        if not 'rtp':
            raise RuntimeError("Failed to locate .rtp files from a local Gromacs installation, please specify"
                               " a custom .rtp file through rtp='/path/to/rtp/file'")
        orig = self.get_atom('resid {} and name CA'.format(resid))
        mutant = gml.ProteinMutant(orig.resname, target)
        targ = mutant.target_3l
        self.top.print("\n  Mutating residue {} (resid {}) into {}\n".format(orig.resname, resid, targ))
        atoms_add, hooks, _, _, extra_bonds, afters = mutant.atoms_to_add()
        atoms_remove = mutant.atoms_to_remove()
        types, charges, dihedrals, impropers, improper_type = self.parse_rtp(rtp, remember=True)
        # some residue-specific modifications here
        if targ == 'HIS':
            targ = 'HSD' if ('HSD', 'CA') in types.keys() else 'HID'
        elif targ == 'GLY':
            self.get_atom('resid {} and name HA'.format(resid)).atomname = 'HA1'
        if orig.resname == 'GLY':
            self.get_atom('resid {} and name HA1'.format(resid)).atomname = 'HA'
        impropers_to_add = []
        impr_sub = self.get_subsection('dihedrals')
        atoms_sub = self.get_subsection('atoms')
        # first remove all unwanted atoms
        for at in atoms_remove:
            equivalents = {'OG': 'OG1', 'HG': 'HG1', 'HG1': 'HG11', 'HG2': 'HG12', 'HG3': 'HG13', 'CG': 'CG1',
                           'CD': 'CD1', 'HD': 'HD1', 'HD1': ['HD11', 'HE2'], 'HD2': 'HD12', 'HD3': 'HD13'}
            self.top.print("Removing atom {} from resid {} in topology".format(at, resid))
            atnum = None
            try:
                atnum = self.get_atom('resid {} and name {}'.format(resid, at)).num
            except RuntimeError:
                equivs = equivalents[at] if isinstance(equivalents[at], list) else [equivalents[at]]
                for equiv in equivs:
                    try:
                        atnum = self.get_atom('resid {} and name {}'.format(resid, equiv)).num
                    except RuntimeError:
                        continue
                    else:
                        break
                else:
                    if atnum is None:
                        raise RuntimeError(f"Couldn't find any of the following: {at}, {', '.join(equivs)}")
            self.del_atom(atnum, del_in_pdb=False)
        for atom_add, hook, aft in zip(atoms_add, hooks, afters):
            self.top.print("Adding atom {} to resid {} in topology".format(atom_add, resid))
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
            xnum = self.get_atom(xsel).num
            ynum = self.get_atom(ysel).num
            self.add_bond(xnum, ynum)
        atoms_sub.get_dicts(force_update=True)
        # looking for new impropers
        for imp in impropers_to_add:
            if set(imp).intersection(set(atoms_add)):
                numbers = []
                for at in imp:
                    if at.startswith('-'):
                        rid = resid - 1
                        atx = at[1:]
                    elif at.startswith('+'):
                        rid = resid + 1
                        atx = at[1:]
                    else:
                        rid = resid
                        atx = at
                    numbers.append(atoms_sub.name_to_num[(rid, atx)])
                # TODO if improper has extra params, add them here
                new_str = '{:5d} {:5d} {:5d} {:5d} {:>5s}\n'.format(*numbers, improper_type)
                impr_sub.add_entry(gml.EntryBonded(new_str, impr_sub),
                                   position=1+[n for n, e in enumerate(impr_sub) if isinstance(e, gml.EntryBonded)][-1])
        # repeating the mutation in the structure
        if mutate_in_pdb and self.top.pdb:
            pdb_atoms = self._match_pdb_to_top(self.get_atom('resid {} and name CA'.format(resid)).num)
            pdb_chains = [self.top.pdb.atoms[at].chain for at in pdb_atoms]
            if len(pdb_atoms) == 1:
                chain = '' if pdb_chains[0] == ' ' else pdb_chains[0]
                self.top.pdb.mutate_protein_residue(resid, target, chain)
            elif len(pdb_atoms) > 1:
                if any([pdb_chains[0] == pdb_chains[i] for i in range(1, len(pdb_chains))]):
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

    def parse_rtp(self, rtp: str, remember: bool = False) -> (dict, dict, dict, dict, dict):
        """
        Reads an .rtp file to extract molecule definitions, separating them into
        dictionaries for: types, charges, impropers, dihedrals, bondedtypes
        :param rtp: str, path to the .rtp file
        :param remember: bool, whether to reuse .rtps selected before
        :return: tuple of dict, each containing atom name : relevant parameter matching
        (ordered: types, charges, dihedrals, impropers, bondedtypes)
        """
        # TODO check against amber/ILDN
        if self.top.rtp and not rtp and remember:
            print("Using previously selected .rtp file")
            return self.top.rtp['typedict'], self.top.rtp['chargedict'], self.top.rtp['dihedrals'], \
                   self.top.rtp['impropers'], self.top.rtp['bondedtypes']
        chargedict, typedict = {}, {}
        impropers, dihedrals = {}, {}
        bondedtypes = 0
        rtp_cont = [line for line in open(rtp) if not line.strip().startswith(';')]
        resname = None
        reading_atoms = False
        reading_impropers = False
        reading_dihedrals = False
        reading_bondedtypes = False
        for line in rtp_cont:
            if line.strip().startswith('[') and line.strip().split()[1] not in ['bondedtypes', 'atoms', 'bonds',
                                                                                'dihedrals', 'impropers']:
                resname = line.strip().split()[1]
            if line.strip().startswith('[') and line.strip().split()[1] == 'atoms':
                reading_atoms = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'atoms':
                reading_atoms = False
            if line.strip().startswith('[') and line.strip().split()[1] == 'dihedrals':
                reading_dihedrals = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'dihedrals':
                reading_dihedrals = False
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
            if len(line.strip().split()) > 3 and resname is not None and reading_dihedrals:
                if resname not in dihedrals.keys():
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
        self.top.rtp['dihedrals'] = dihedrals
        self.top.rtp['impropers'] = impropers
        # TODO we need a method for explicitly setting just a single dihedral (should be easy)
        self.top.rtp['bondedtypes'] = bondedtypes
        return self.top.rtp['typedict'], self.top.rtp['chargedict'], self.top.rtp['dihedrals'], \
            self.top.rtp['impropers'], self.top.rtp['bondedtypes']

    def update_dicts(self):
        self.get_subsection('atoms').get_dicts(force_update=True)

    def list_atoms(self):
        for atom in self.atoms:
            print(str(atom).strip())

    def list_bonds(self, by_types: bool = False, by_params: bool = False, returning: bool = False):
        return self._list_bonded('bonds', by_types, by_params, returning)

    def list_angles(self, by_types: bool = False, by_params: bool = False, returning: bool = False):
        return self._list_bonded('angles', by_types, by_params, returning)

    def list_dihedrals(self, by_types: bool = False, by_params: bool = False, returning: bool = False,
                       interaction_type: Optional[list]=None):
        return self._list_bonded('dihedrals', by_types, by_params, returning, interaction_type)

    def to_rtp(self, outname='out.rtp'):
        if os.path.exists(outname):
            mode = 'a'
            print('Found an existing .rtp file with requested name, will append to it')
        else:
            mode = 'w'
        if not all([self.atoms[0].resname == a.resname for a in self.atoms]):
            raise RuntimeError(f'Not all atoms in molecule {self.mol_name} have the same residue name, '
                               f'fix this to write an .rtp entry')
        resname = self.atoms[0].resname
        atoms = [(a.atomname, a.type, a.charge, n) for n, a in enumerate(self.atoms, 1)]
        bonds = self.list_bonds(returning=True)
        impropers = self.list_dihedrals(returning=True, interaction_type=["2", "4"])
        with open(outname, mode) as outfile:
            outfile.write(f"[ {resname} ]\n")
            outfile.write(f"  [ atoms ]\n")
            for at in atoms:
                outfile.write(f"  {at[0]:6s} {at[1]:6s} {at[2]:8.4f} {at[3]:5d}\n")
            outfile.write(f"  [ bonds ]\n")
            for bd in bonds:
                outfile.write(f"  {bd[0]:6s} {bd[1]:6s}\n")
            if impropers:
                outfile.write(f"  [ impropers ]\n")
                for imp in impropers:
                    outfile.write(f"  {imp[0]:6s} {imp[1]:6s} {imp[2]:6s} {imp[3]:6s}\n")

    def _list_bonded(self, term: str, by_types: bool = False, by_params: bool = False, returning: bool = False,
                     interaction_type: Optional[list] = None):
        self.update_dicts()
        subsection = self.get_subsection(term)
        tried_adding = False
        returnable = []
        formatstring = {'bonds': "{:>5s} {:>5s}", 'angles': "{:>5s} {:>5s} {:>5s}",
                        'dihedrals': '{:>5s} {:>5s} {:>5s} {:>5s}'}
        for entry in subsection:
            if isinstance(entry, gml.EntryBonded):
                entry.read_types()
                if not by_params:
                    extra = ''
                    params = []
                else:
                    if not entry.params_state_a and not tried_adding:
                        self.top.add_ff_params()
                        tried_adding = True
                    extra = '{:>12.5f} ' * len(entry.params_state_a)
                    params = entry.params_state_a
                if interaction_type is not None:
                    if entry.interaction_type not in interaction_type:
                        continue
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

    def alch_h_to_ch3(self, resid: int, orig_name: str, basename: str, ctype: Optional[str] = None,
                      htype: Optional[str] = None, ccharge: Optional[float] = None, hcharge: float = 0.09,
                      dummy_type: str = 'DH', add_in_pdb: bool = True):
        """
        A generic routine to change a hydrogen atom into a methyl group
        :param resid: int, ID of the residue in which the change is to be made
        :param orig_name: str, name of the hydrogen to be modified
        :param basename: str, an atom name for the extra atoms that will be modified with indices
        :param ctype: str, atomtype for the methyl carbon
        :param htype: str, atomtype for the methyl hydrogens
        :param ccharge: float, charge for the methyl carbon
        :param hcharge: float, charge for the methyl hydrogens
        :param dummy_type: str, atomtype for the dummies
        :param add_in_pdb: bool, whether to modify the associated Pdb object
        :return: None
        """
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
        orig = self.get_atom('resid {}'.format(resid))
        ccharge = ccharge if ccharge is not None else round(orig.charge - 0.27, 4)
        atoms_add, hooks = [basename.replace('C', 'H') + str(i) for i in range(3)], 3 * [orig_name]
        for n, atom_add_hook in enumerate(zip(atoms_add, hooks), 1):
            atom_add, hook = atom_add_hook
            self.top.print("Adding atom {} to resid {} in the topology".format(atom_add, resid))
            hooksel = 'resid {} and name {}'.format(resid, orig_name)
            hnum = self.get_atom(hooksel).num
            atnum = hnum + n
            self.add_atom(atnum, atom_add, atom_type=dummy_type, charge=0, resid=resid, resname=orig.resname, mass=1.008)
            self.add_bond(hnum, atnum)
            self.gen_state_b(atomname=atom_add, resid=resid, new_type=htype, new_charge=hcharge, new_mass=1.008)
        self.gen_state_b(atomname=orig_name, resid=resid, new_type=ctype, new_charge=ccharge, new_mass=12.0)
        if add_in_pdb and self.top.pdb:
            if len(self.top.system) > 1 or self.top.system[0][1] > 1:
                raise RuntimeError("Adding groups in PDB only supported for systems containing one molecule")
            bonds = self.list_bonds(returning=True)
            hook = [j for i in bonds for j in i if orig_name in i and orig_name != j][0]
            aligns = [j for i in bonds for j in i if hook in i and hook != j and orig_name != j]
            aftnr = self.select_atom('resid {} and name {}'.format(resid, orig_name))
            for n, aliat in enumerate(zip(aligns, atoms_add), 1):
                ali, at = aliat
                self.top.pdb.insert_atom(aftnr+n, self.top.pdb.atoms[aftnr],
                                         atomsel='resid {} and name {}'.format(resid, at),
                                         hooksel='resid {} and name {}'.format(resid, orig_name), bondlength=1.1,
                                         p1_sel='resid {} and name {}'.format(resid, ali),
                                         p2_sel='resid {} and name {}'.format(resid, hook), atomname=at)

    def add_dummy_def(self, dummy_type: str):
        params = self.top.parameters
        atomtypes = params.get_subsection('atomtypes')
        dummy_entries = [e for e in atomtypes if isinstance(e, gml.EntryParam) and e.types[0] == dummy_type]
        if not dummy_entries:
            atomtypes.add_entry(gml.EntryParam('   {}     0        1.008  0.0000  A  0.000000000000  0.0000  '
                                               '\n'.format(dummy_type), atomtypes))

    def add_vs2(self, atom1: int, atom2: int, fraction: float, atomname: str, atomtype: str):
        """
        Adds an interpolation-type virtual site (between two atoms) at the end of the molecule.
        If the respective section does not exist in the topology, it will be created.
        :param atom1: int, serial numbers of the first atom
        :param atom2: int, serial numbers of the second atom
        :param fraction: float, at which point between the 1st (0) and 2nd (1) atom the VS will be created
        :param atomname: str, name of the new atom
        :param atomtype: str, type of the new atom
        :return: None
        """
        try:
            ssect = self.get_subsection('virtual_sites2')
        except KeyError:
            self.subsections.append(self._yield_sub(['[ virtual_sites2 ]']))
            ssect = self.get_subsection('virtual_sites2')
        ref = self.atoms[atom1-1]
        self.add_atom(len(self.atoms) + 1, atomname, atomtype, resid=ref.resid, resname=ref.resname, mass=0.0)
        ssect.add_entry(gml.EntryBonded(f'{atom1} {atom2} {len(self.atoms)} 1 {fraction}\n', ssect))

    def add_constraint(self, atom1: int, atom2: int, distance: float):
        """
        Adds a constraint
        :param atom1: int, serial numbers of the first atom
        :param atom2: int, serial numbers of the second atom
        :param distance: float, the distance (in nm) at which the pair will be constrained
        :return: None
        """
        try:
            ssect = self.get_subsection('constraints')
        except KeyError:
            self.subsections.append(self._yield_sub(['[ constraints ]']))
            ssect = self.get_subsection('constraints')
        ssect.add_entry(gml.EntryBonded(f'{atom1} {atom2} 2 {distance}\n', ssect))

    def make_stateB_dummy(self, resid: int, orig_name: str, dummy_type: str = 'DH'):
        self.add_dummy_def(dummy_type)
        self.gen_state_b(atomname=orig_name, resid=resid, new_type=dummy_type, new_charge=0, new_mass=1.008)

    def solute_tempering(self, temperatures: list, selection: str = None):
        """
        Prepares a modified topology for REST2 simulations
        and writes the respective .top files
        :param temperatures: list of float, solute temperatures for replica exchange
        :param selection: str or None, can specify the selection to which the modification should be restricted
        :return: None
        """
        self.top.explicit_defines()
        for n, t in enumerate(temperatures):
            self.top.print(f'generating topology for effective temperature of {t} K...')
            mod = deepcopy(self.top)
            mod.get_molecule(self.mol_name).scale_rest2_charges(temperatures[0]/t, selection)
            mod.get_molecule(self.mol_name).scale_rest2_bonded(temperatures[0]/t, selection)
            mod.save_top(self.top.fname.replace('.top', f'-rest{temperatures[0]/t:.3f}.top'))

    def scale_rest2_bonded(self, gamma: float, selection: str = None):
        """
        Modifies bonded terms for the REST2 protocol
        :param gamma: float, the scaling parameter from REST2
        :param selection: str or None, can specify the selection to which the scaling should be restricted
        :return: None
        """
        # get a list of atomtypes to clone
        sel = self.get_atoms(selection) if selection is not None else self.atoms
        types = {at.type for at in sel}
        for tp in types:
            self.top.parameters.clone_type(tp)
        for sub in ['atomtypes', 'bondtypes', 'angletypes', 'dihedraltypes', 'nonbond_params', 'pairtypes']:
            try:
                for subsect in self.top.parameters.get_subsections(sub):
                    for ent in subsect.entries_param:
                        if any([tp.startswith('y') for tp in ent.types]):
                            ent.params[1] *= gamma
            except KeyError:
                pass
        try:
            for ent in self.top.parameters.get_subsection('cmaptypes').entries_param:
                if any([tp.startswith('y') for tp in ent.types]):
                    ent.params = [float(x)*gamma for x in ent.params]
        except KeyError:
            pass

    def scale_rest2_charges(self, gamma: float, selection: str = None):
        """
        Modifies charges for the REST2 protocol
        :param gamma: float, the scaling parameter from REST2
        :param selection: str or None, can specify the selection to which the scaling should be restricted
        :return: None
        """
        sel = self.get_atoms(selection) if selection is not None else self.atoms
        # scaling charges
        for a in self.atoms:
            if a in sel:
                a.charge = round(a.charge * gamma**0.5, 4)
                a.type = 'y' + a.type

    def alchemical_proton(self, resid: int, rtp: Optional[str] = None, b_is_protonated: bool = False):
        """
        Creates an alchemical residue (starting from ASP or GLU) where the B-state
        corresponds to a protonated variant of that residue
        :param resid: int, number of the residue to modify
        :param rtp: str, path to the aminoacids.rtp or merged.rtp file (optional)
        :param b_is_protonated: bool, whether to make the B-state protonated; default is False (state A is protonated)
        :return: None
        """
        atoms = self.get_atoms(f'resid {resid}')
        resname = atoms[0].resname
        if resname not in ['ASP', 'GLU', 'GLUP', 'ASPP']:
            raise RuntimeError("So far only available for residues ASP/GLU, and ASPP/GLUP")
        mut_dict = {'ASP': 'B', 'GLU': 'J'}
        deprot_dict = {'GLUP': 'GLU'}
        if resname in ['ASP', 'GLU']:
            self.mutate_protein_residue(resid, mut_dict[resname])
        else:
            for atom in self.get_atoms(f'resid {resid}'):
                atom.resname = deprot_dict[resname]
            resname = deprot_dict[resname]
        atoms = self.get_atoms(f'resid {resid}')
        rtp = self.find_rtp(rtp)
        types, charges, dihedrals, impropers, improper_type = self.parse_rtp(rtp)
        for atom in atoms:
            try:
                atom.type_b = types[(resname, atom.atomname)]
                atom.charge_b = charges[(resname, atom.atomname)]
                atom.mass_b = atom.mass
            except KeyError:
                atom.type_b = 'DH'
                atom.charge_b = 0.0
                atom.mass_b = 1.008
        if 'DH' not in self.top.defined_atomtypes:
            subsect = self.top.parameters.get_subsection('atomtypes')
            subsect.add_entry(gml.EntryParam('DH  0  0.0  0.0  A   0.0  0.0', subsect))
        if b_is_protonated:
            self.swap_states(resid=resid)
        self.update_dicts()
        self.add_ff_params()

    def recalc_qtot(self):
        """
        Puts the qtot (cumulative charge of the molecule) in the comment
        of an atom entry in the topology (done by pdb2gmx by default)
        :return: None
        """
        charge = 0
        for atom in self.atoms:
            charge += atom.charge
            atom.comment = f'qtot {charge:.3f}'


class SectionParam(Section):
    """
    This class should wrap together sections such as [ bondtypes ],
    [ atomtypes ], [ pairtypes ] etc. and have methods designed to
    facilitate the search of matching params
    """
    
    def __init__(self, content_list: list, top: "gml.Top"):
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
        subsection_labels = [sub.header for sub in self.subsections]
        duplicated_subsections = list({label for label in subsection_labels if subsection_labels.count(label) > 1})
        for sub in duplicated_subsections:
            subsections_to_merge = [s for s in self.subsections if s.header == sub and not s.conditional]
            if not subsections_to_merge:
                continue
            self.top.print(f"Merging sections {subsections_to_merge} together")
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
                sub.sort()  # TODO if two have same periodicity & atoms, remove one with 0-s (PMX)

    def add_dummy_def(self, dummy_type: str):
        atomtypes = self.get_subsection('atomtypes')
        dummy_entries = [e for e in atomtypes if isinstance(e, gml.EntryParam) and e.types[0] == dummy_type]
        if not dummy_entries:
            atomtypes.add_entry(gml.EntryParam('   {}     0        1.008  0.0000  A  0.000000000000  0.0000  '
                                               '\n'.format(dummy_type), atomtypes))

    def find_used_ff_params(self, section: str = 'all') -> list:
        """
        Finds FF parameters that are used by the system
        :param section: str, for which section should the search be performed ('all' or 'atomtypes', '...')
        :return: list of EntryParams that are used by the system's topology
        """
        used_params = []
        if section == 'all':
            subsections_to_add = ['atomtypes', 'pairtypes', 'nonbond_params', 'constrainttypes']
            # TODO check what with other sections?
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

    def fix_zero_periodicity(self):
        """
        Recent versions of OpenMM require non-zero dihedral periodicities,
        so this function changes 0s to 1s (assuming the force constant is 0 too)
        :return: None
        """
        subsections = self.get_subsections('dihedraltypes')
        for sub in subsections:
            for entry in sub.entries_param:
                if entry.params[-1] == 0 and entry.interaction_type in '149':
                    assert entry.params[-2] == 0
                    entry.params[-1] = 1

    def clone_type(self, atomtype: str, prefix: str = 'y', new_type: Optional[str] = None):
        """
        Generates an exact type of a selected atomtype,
        preserving all interactions with other types
        :param atomtype: str, atomtype to be duplicated
        :param prefix: str, new name will be generated as prefix + original atomtype
        :param new_type: str, directly specify the new name (optional)
        :return: None
        """
        new_atomtype = prefix + atomtype if new_type is None else new_type
        existing_types = {e.types[0] for e in self.top.parameters.get_subsection('atomtypes').entries_param}
        if new_atomtype in existing_types:
            raise RuntimeError(f"Type {new_atomtype} already exists")
        for sub in self.subsections:
            to_add = []
            for ent in sub:
                if isinstance(ent, gml.EntryParam) and atomtype in ent.types and new_atomtype not in ent.types:
                    to_add.append(ent)
            for entry in to_add:
                newlines = self.gen_clones(entry, atomtype, new_atomtype)
                sub.add_entries([gml.EntryParam(line, sub) for line in newlines])
        self.sort_dihedrals()
        self._remove_symm_dupl(new_atomtype)

    def clean_unused(self, used_params: str, section: str = 'all'):
        """
        Cleans up FF parameters that are not used by the given system
        :param used_params: list of str, identifiers of parameters we already know are being used
        :param section: str, which sections to clean up ('all' for all of them or 'bonds', 'angles', 'dihedrals', etc.)
        :return: None
        """
        matchings = {'bonds': 'bondtypes', 'angles': 'angletypes', 'dihedrals': 'dihedraltypes',
                     'atomtypes': 'atomtypes', 'pairtypes': 'pairtypes', 'nonbond_params': 'nonbond_params',
                     'constrainttypes': 'constrainttypes', 'cmap': 'cmaptypes',
                     'implicit_genborn_params': 'implicit_genborn_params'}
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
        atomtypes_used = {e.type for mol in self.top.molecules for e in mol.atoms
                          if isinstance(e, gml.EntryAtom)}
        atomtypes_b_used = {e.type_b for mol in self.top.molecules for e in mol.atoms
                            if isinstance(e, gml.EntryAtom) and e.type_b}
        atomtypes_used.union(atomtypes_b_used)
        ssect = self.get_subsection('atomtypes')
        new_entries = []
        for entry in ssect.entries:
            if not isinstance(entry, gml.EntryParam) or entry.types[0] in atomtypes_used:
                new_entries.append(entry)
        ssect.entries = new_entries

    def _remove_symm_dupl(self, new_atomtype: str):
        for sub in self.subsections:
            if 'dihedral' in sub.header:
                sub._remove_symm(new_atomtype)

    def get_opt_dih(self, types: bool = False) -> list:
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and '9' in sub.prmtypes][0]
        return ss.get_opt_dih(types)

    def get_opt_dih_indices(self) -> list:
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and '9' in sub.prmtypes][0]
        return ss.get_opt_dih_indices()

    def set_opt_dih(self, values):
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and '9' in sub.prmtypes][0]
        ss.set_opt_dih(values)

    def add_nbfix(self, type1: str, type2: str, mod_sigma: float = 0.0, mod_epsilon: float = 0.0,
                  action_default: str = 'x', new_sigma: Optional[float] = None, new_epsilon: Optional[float] = None):
        """
        Generates NBFIX entries for the chosen pair of atomtypes by modifying the current
        Lorentz-Berthelot rules, or an existing NBFIX
        :param type1: str, name of the 1st atomtype in the pair
        :param type2: str, name of the 2nd atomtype in the pair
        :param mod_sigma: float, by how much to increase the LJ sigma (in nm)
        :param mod_epsilon: float, by how much to increase the LJ epsilon (in kJ/mol)
        :param new_sigma: float, if specified instead of mod_sigma sets the value directly (in nm)
        :param new_epsilon: float, if specified instead of mod_epsilon sets the value directly (in kJ/mol)
        :param action_default: str, what to do if an NBFIX already exists (check prompt if that happens)
        :return: None
        """
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
        new_sig = 0.5 * (sigma1 + sigma2) + mod_sigma if new_sigma is None else new_sigma
        new_eps = (eps1 * eps2) ** 0.5 + mod_epsilon if new_epsilon is None else new_epsilon
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
                        new_sig = entry.params[0] + mod_sigma if new_sigma is None else new_sigma
                        new_eps = entry.params[1] + mod_epsilon if new_epsilon is None else new_epsilon
                        comment = entry.comment
                    nbsub.remove_entry(entry)
        entry_line = "{} {} 1 {} {} ; sigma chg by {}, eps chg by {} {}".format(type1, type2, new_sig, new_eps,
                                                                                mod_sigma, mod_epsilon, comment)
        nbsub.add_entry(gml.Subsection.yield_entry(nbsub, entry_line))

    def edit_atomtype(self, atomtype: str, mod_sigma: float = 0.0, mod_epsilon: float = 0.0,
                      new_sigma: Optional[float] = None, new_epsilon: Optional[float] = None):
        """
        Modifies the values of sigma or epsilon for a chosen atomtype
        :param atomtype: str, type to be edited
        :param mod_sigma: float, by how much should sigma be changed (in nm)
        :param mod_epsilon: float, by how much should epsilon be changed (in kJ/mol)
        :param new_sigma: float, if specified instead of mod_sigma sets the value directly (in nm)
        :param new_epsilon: float, if specified instead of mod_epsilon sets the value directly (in kJ/mol)
        :return: None
        """
        atp = self.get_subsection('atomtypes')
        for entry in atp:
            if isinstance(entry, gml.EntryParam):
                if entry.types[0] == atomtype:
                    entry.params[0] = entry.params[0] + mod_sigma if new_sigma is None else new_sigma
                    entry.params[1] = entry.params[1] + mod_epsilon if new_epsilon is None else new_epsilon
                    # TODO adjust comment if using new_[se]
                    entry.comment = f"; sigma chg by {mod_sigma}, eps chg by {mod_epsilon} {entry.comment}"
                    return
        raise RuntimeError(f"Couldn't find type {atomtype}, check your topology")

    def add_atomtype(self, atomtype: str, mass: float, sigma: float, epsilon: float, action_default: str = 'x',
                     atomic_number: Optional[int] = None):
        """
        Adds a new atomtype, defining mass, sigma, epsilon, and (optionally) atomic number
        :param atomtype: str, name of the type
        :param mass: float, mass of the atom
        :param sigma: float, sigma of the Lennard-Jones potential (in nm)
        :param epsilon: float, epsilon of the Lennard-Jones potential (in kJ/mol)
        :param action_default: str, what to do if an entry already exists (check prompt if that happens)
        :param atomic_number: int, optional
        :return: None
        """
        atnum = atomic_number if atomic_number is not None else 0
        atp = self.get_subsection('atomtypes')
        if atomtype in self.top.defined_atomtypes:
            action = action_default
            while action not in 'rt':
                action = input("An entry already exists, shall we replace it (r) or terminate (t)?")
            if action == 't':
                return
            else:
                entry = [ent for ent in atp.entries_param if ent.types[0] == atomtype][0]
                print(f"Overwriting entry for atomtype {atomtype} with mass {mass}, sigma {sigma}, epsilon {epsilon}")
                del entry
        atp.add_entry(gml.EntryParam(f'{atomtype} {atnum} {mass} 0.0000 A {sigma} {epsilon}', subsection=atp))

    def add_bonded_param(self, types: tuple, params: list, interaction_type: int, action_default: str = 'x'):
        """
        # TODO fill here
        :param types:
        :param params:
        :param interaction_type:
        :param action_default:
        :return:
        """
        subsection_dict = {(2, 1): 'bondtypes', (3, 1): 'angletypes', (4, 1): 'dihedraltypes', (4, 9): 'dihedraltypes',
                           (4, 4): 'dihedraltypes'}
        subs = self.get_subsection(subsection_dict[(len(types), interaction_type)])
        matching = [ent for ent in subs.entries_param if ent.types == types or ent.types == types[::-1]]
        if matching:
            action = action_default
            while action not in 'rt':
                action = input("An entry already exists, shall we replace it (r) or terminate (t)?")
            if action == 't':
                return
            else:
                while matching:
                    del matching[-1]
        if len(params) <= 3:
            subs.add_entry(gml.EntryParam(f'{" ".join(types)} {str(interaction_type)} '
                                          f'{" ".join([str(x) for x in params])}', subsection=subs))
        elif interaction_type == 9:
            assert len(params) % 3 == 0
            for i in range(len(params) // 3):
                param_slice = params[i*3: (i+1)*3]
                subs.add_entry(gml.EntryParam(f'{" ".join(types)} {str(interaction_type)} '
                                              f'{" ".join([str(x) for x in param_slice])}', subsection=subs))

    @staticmethod
    def gen_clones(entry: "gml.EntryParam", atomtype: str, new_atomtype: str) -> list:
        """
        Copies entry that contains an atomtype to be cloned, taking into
        account possible multiple occurences of that type and generating
        all possible combinations/permutations of the changes
        :param entry: an EntryParam instance, the entry to be copied
        :param atomtype: str, the atomtype that is being cloned
        :param new_atomtype: str, the new atomtype being created
        :return: list of str, new lines to be added
        """
        lines = []
        nchanges = entry.types.count(atomtype)
        changes = []
        for i in range(nchanges):
            changes.extend(SectionParam.gen_combs(nchanges, i + 1))
        for mod in changes:
            lines.append(SectionParam.mod_types(entry, mod, new_atomtype, atomtype))
        return lines

    @staticmethod
    def gen_combs(count: int, tuples: int) -> list:
        return list(combinations(range(count), tuples))

    @staticmethod
    def mod_types(entry: "gml.EntryParam", mods: list, new_atomtype: str, atomtype: str) -> str:
        """
        Modifies an EntryParam to include the new atomtype instead of the old one
        in places specified by the directive in mods
        :param entry: an EntryParam instance, the entry to be copied
        :param mods: list of lists, specifies where to put the modification
        :param new_atomtype: str, the new atomtype being created
        :param atomtype: str, the atomtype that is being cloned
        :return: str, the cloned line
        """
        line = str(entry)
        lentype = len(atomtype)
        for num in mods[::-1]:
            indices = [i for i in range(len(line) - len(atomtype) + 1)
                       if line[i:i + len(atomtype)] == atomtype and (i == 0 or line[i-1].isspace())
                       and (i+len(atomtype) == len(line) or line[i+len(atomtype)].isspace())]
            line = line[:indices[num]] + new_atomtype + line[indices[num]+lentype:]
        return line
