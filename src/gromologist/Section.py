"""
Module: Section.py
Author: Miłosz Wieczór <milosz.wieczor@irbbarcelona.org>
License: GPL 3.0

Description:
    This module represents a portion of the topology corresponding to one of the
    setup directives, a single molecule, or the force field parameter set.

Contents:
    Classes:
        Section:
            Base class, represents any section and implements basic
            bookkeeping + getters/setters
        SectionMol:
            Subclass, represents a molecule and allows for many molecule-level
            operations, from editing atoms to making molecules alchemical
        SectionParam:
            Subclass, represents the entire parameter set associated with the system;
            only one is allowed per topology, allows for general FF adjustments

Usage:
    This module is intended to be imported as part of the library. It is not
    meant to be run as a standalone script. Example:

        import gromologist as gml
        t = gml.Top("complex.top")
        print(t.molecules)  # molecules is an alias for a list of SectionMol's

Notes:
    All operations on molecules assume 1-based indexing, consistent with Gromacs'
    topology numbering for atoms
"""


import os
from time import time
from itertools import product, combinations
from typing import Optional, Union
from functools import reduce
from copy import deepcopy
from glob import glob

import gromologist as gml
from gromologist import Subsection


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
                    self.top.print(f"Subection {str(self.subsections[-1])} recognized as conditional, "
                                   f"will not be merged")
                elif isinstance(self.subsections[-1], gml.SubsectionBonded) \
                        or isinstance(self.subsections[-1], gml.SubsectionAtom):
                    self.top.print(f"Subection {str(self.subsections[-1])} in molecule "
                                   f"{self.subsections[-1].section.subsections[0].molname} recognized as conditional")
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
                      'virtual_sitesn', 'virtual_sites2', 'virtual_sites3', 'constraints', 'pairs_nb'}:
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

    def init_subsection(self, subsection_name) -> None:
        self.subsections.append(self._yield_sub([f'[ {subsection_name} ]\n']))

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
        """
        Returns gml.EntryAtom entries corresponding to all atoms in the molecule
        :return: list of gml.EntryAtom
        """
        return self.get_subsection('atoms').entries_atom

    @property
    def natoms(self) -> int:
        """
        Returns the number of atoms in the molecule
        :return: int, number of atoms
        """
        return len(self.atoms)

    @property
    def nmols(self) -> int:
        """
        Returns the count of the given molecule in the complete system, according to the [ molecules ] section
        :return: int, number of molecule copies
        """
        return sum([mol_count[1] for mol_count in self.top.system if mol_count[0] == self.mol_name])

    @property
    def charge(self) -> float:
        """
        Calculates the total charge of the molecule
        :return: float, total charge
        """
        return sum([a.charge for a in self.atoms])

    @property
    def mass(self) -> float:
        """
        Calculates the total mass of the molecule
        :return: float, total mass
        """
        return sum([a.mass for a in self.atoms])

    @property
    def bonds_section(self) -> "gml.SubsectionBonded":
        """
        Returns the [ bonds ] section of the topology
        :return: gml.SubsectionBonded
        """
        return self.get_subsection('bonds')

    @property
    def angles_section(self) -> "gml.SubsectionBonded":
        """
        Returns the [ angles ] section of the topology
        :return: gml.SubsectionBonded
        """
        return self.get_subsection('angles')

    @property
    def dihedrals_section(self) -> "gml.SubsectionBonded":
        """
        Returns the [ dihedrals ] section of the topology
        :return: gml.SubsectionBonded
        """
        return self.get_subsection('dihedrals')

    @property
    def residues(self) -> list:
        """
        Returns the list of residues defined in the topology (labeled as resname-resid pairs)
        :return: list of str
        """
        resid = None
        reslist = []
        for at in self.atoms:
            if at.resid != resid:
                resid = at.resid
                reslist.append(f'{at.resname}-{at.resid}')
        return reslist

    @property
    def atomtypes(self) -> list:
        """
        Returns a list of all the atom types in the molecule
        :return: list of str, all atom types used in the molecule
        """
        allt = {atom.type for atom in self.atoms}.union({atom.type_b for atom in self.atoms if atom.type_b is not None})
        return sorted(list(allt))

    def _merge(self) -> None:
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

    def _remove_subsection(self, header) -> None:
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

    def set_type(self, type_to_set: Optional[str] = None, atomname: Optional[str] = None, resname: Optional[str] = None,
                 resid: Optional[int] = None, atomtype: Optional[str] = None, prefix: Optional[str] = None) -> None:
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
        if type_to_set is None and prefix is None:
            raise RuntimeError("Set either type_to_set or prefix")
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
        if atomtype is None:
            types = list({a.type for a in self.atoms})
        else:
            if isinstance(atomtype, list) or isinstance(atomtype, tuple):
                types = atomtype
            else:
                types = [atomtype]
        if atomname is None:
            names = list({a.atomname for a in self.atoms})
        else:
            if isinstance(atomname, list) or isinstance(atomname, tuple):
                names = atomname
            else:
                names = [atomname]
        resnames, names, resids, types = set(resnames), set(names), set(resids), set(types)
        for a in self.atoms:
            if a.resname in resnames and a.atomname in names and a.resid in resids and a.type in types:
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
        """
        Key getter for atoms, returns a list of gml.EntryAtom instances
        when queried with a selection
        :param selection_string: str, selection/query
        :return: list of gml.EntryAtom, all entries that correspond to the requested atoms
        """
        return [self.atoms[i] for i in self.select_atoms(selection_string)]

    def get_atom(self, selection_string: str) -> "gml.EntryAtom":
        """
        Key getter for atoms, returns a single gml.EntryAtom instance
        when queried with a selection; the selection has to point to a unique atom
        :param selection_string: str, selection/query
        :return: gml.EntryAtom, an entry corresponding to the requested atom
        """
        return self.atoms[self.select_atom(selection_string)]

    def print_molecule(self):
        """
        Prints a list of atom in the molecule, useful when dynamically editing the molecule
        :return: None
        """
        sub = self.get_subsection('atoms')
        for entry in sub:
            print(str(entry), end='')

    @property
    def is_alchemical(self) -> bool:
        """
        Whether the molecule contains any alchemical atoms/residues
        :return: bool
        """
        sect = self.get_subsection('atoms')
        for ent in sect.entries_atom:
            if ent.type_b is not None:
                return True
        return False

    @property
    def is_water(self) -> bool:
        """
        Whether the molecule is a water residue
        :return: bool
        """
        return bool(self.get_atoms('water'))

    @property
    def is_protein(self) -> bool:
        """
        Whether the molecule is a protein
        :return: bool
        """
        return bool(self.get_atoms('protein'))

    @property
    def is_nucleic(self) -> bool:
        """
        Whether the molecule is a nucleic acid
        :return: bool
        """
        return bool(self.get_atoms('nucleic'))

    def _patch_alch(self) -> None:
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
                self._nullify_bonded(*atoms[2*pair:2*pair+2], 'angles')
                self._nullify_bonded(*atoms[2 * pair:2 * pair + 2], 'dihedrals')
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

    def _nullify_bonded(self, atom1, atom2, subsection, stateB=True) -> None:
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

    def offset_numbering(self, offset: int, startfrom: int = 0) -> None:
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

    def offset_residues(self, offset: int, startfrom: int = 0) -> None:
        """
        Offsets residue numbering starting from a specified position;
        useful when trying to match a system to a reference one
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        for a in self.atoms:
            if a.num >= startfrom:
                a.resid += offset

    def _offset_atoms(self, offset: int, startfrom: int) -> None:
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

    def _offset_params(self, offset: int, startfrom: int) -> None:
        """
        Offsets atomic numbering in all parameter sections,
        e.g., [ bonds ]
        :param offset: int, by how much we wish to offset the numbering
        :param startfrom: int, starting point of the offset
        :return: None
        """
        for sub_name in [s.header for s in self.subsections if s.header != 'atoms']:
            subsections = self.get_subsections(sub_name)
            for subsection in subsections:
                try:
                    for entry in subsection.entries_bonded:
                        entry.atom_numbers = tuple(n + (offset * (n >= startfrom)) for n in entry.atom_numbers)
                except AttributeError:
                    continue

    def add_disulfide(self, resid1: int, resid2: int, other: Optional["gml.SectionMol"] = None,
                      rtp: Optional[str] = None) -> None:
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
        # TODO remove in PDB
        other = self if other is None else other
        s1 = self.get_atom(f'resid {resid1} and name SG')
        h1 = self.get_atom(f'resid {resid1} and name HG HG1')
        s2 = other.get_atom(f'resid {resid2} and name SG')
        h2 = other.get_atom(f'resid {resid2} and name HG HG1')
        assert s1.resname == s2.resname == 'CYS'
        self.del_atom(h1.num)
        other.del_atom(h2.num)
        found = self.find_rtp(rtp)
        types, charges, dihedrals, impropers, improper_type, _ = self.parse_rtp(found)
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
        self._check_correct()

    def add_coordinated_ion(self, resid1: int, resid2: int, other: Optional["gml.SectionMol"] = None,
                            rtp: Optional[str] = None) -> None:
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
        # TODO remove in PDB
        other = self if other is None else other
        c1 = self.get_atom(f'resid {resid1} and name CA')
        if c1.resname not in {'CYS', 'HIS', 'HID', 'HIE', 'HSD', 'HSE'}:
            raise RuntimeError(f"Residue 1 has to be cysteine or histidine, not {c1.resname}")
        x2 = other.get_atom(f'resid {resid2} and name ZN ZN2 FE')
        if x2.resname not in {'ZN', 'ZN2', 'FE', 'HEM', 'HEME'}:
            raise RuntimeError(f"Residue 2 has to be a zinc or iron/heme, not {x2.resname}")
        if c1.resname == 'CYS':
            found = self.find_rtp(rtp)
            types, charges, dihedrals, impropers, improper_type, _ = self.parse_rtp(found)
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
        self._check_correct()

    def convert_to_hydrogens(self, selection: str, neutralize_h: bool = False, target_mol_charge: float = 0.0) -> None:
        """
        A function to fix sliced MM subsystems with broken bonds by converting dangling atoms to hydrogens;
        useful for running comparisons between QM and MM subsystems
        :param selection: str, a selection that defines atoms to be converted to hydrogens
        :param neutralize_h: bool, whether to modify/redistribute the charge to give the subsystem an integer charge
        :param target_mol_charge: float, if neutralizing, this defines what the total charge should be
        :return: None
        """
        to_mod = self.get_atoms(selection)
        if any([a.numbonds > 1 for a in to_mod]):
            bdict = {a.atomname : a.numbonds for a in to_mod}
            raise AttributeError(f'At least one of the selected atoms has more than one binding partner: {bdict}')
        adjustable_hydrogens = []
        for atom in to_mod:
            bound = atom.bound_atoms[0]
            next_bounds = bound.bound_atoms
            bounds_h = [b for b in next_bounds if b.ish]
            assert all([b.mass == bounds_h[0].mass for b in bounds_h]) and \
                   all([b.charge == bounds_h[0].charge for b in bounds_h]) and \
                   all([b.type == bounds_h[0].type for b in bounds_h])
            atom.type = bounds_h[0].type
            atom.charge = bounds_h[0].charge
            atom.mass = bounds_h[0].mass
            self.top.print(f"Setting {atom.atomname} to type {atom.type}, charge {atom.charge}, mass {atom.mass}; "
                           f"renaming to {'Hx' + atom.atomname}")
            atom.atomname = 'Hx' + atom.atomname
            if self.top.pdb:
                matched_h = self._match_pdb_to_top(atom.num)
                matched_heavy = self._match_pdb_to_top(bound.num)
                for mh, mhea in zip(matched_h, matched_heavy):
                    self.top.pdb.reposition_atom_from_hook(atomsel=f'serial {mh}', hooksel=f'serial {mhea}',
                                                           bondlength=1.09, p1_sel=f'serial {mhea}',
                                                           p2_sel=f'serial {mh}')
                    self.top.pdb.get_atom(f'serial {mh}').atomname = atom.atomname[:4]
            for bh in bounds_h + [atom]:
                if bh not in adjustable_hydrogens:
                    adjustable_hydrogens.append(bh)
        if neutralize_h and adjustable_hydrogens:
            excess_charge_per_h = round((self.charge - target_mol_charge) / len(adjustable_hydrogens), 8)
            for adj in adjustable_hydrogens:
                adj.charge -= excess_charge_per_h

    def alchemize(self, resid: int, new_resname: str, rtp: str, matching_atoms: Optional[dict] = None) -> None:
        """
        Creates an alchemical residue from an existing residue, another residue template,
        and a dictionary of corresponding atoms; all non-matching atoms (i.e. different names
        or not included in matching_atoms) will be created/shrunk as dummies
        :param resid: int, number of the residue to be converted
        :param new_resname: str, name of the residue to be merged with
        :param rtp: str, path to the .rtp file that contains the residue template
        :param matching_atoms: dict, matches atoms in the original residue to atoms in the new residue
        :return: None
        """
        types, charges, dihedrals, impropers, bondedtypes, bonds = self.parse_rtp(rtp)
        if new_resname not in {x[0] for x in types.keys()}:
            raise RuntimeError(f"New residue name {new_resname} not found in file {rtp}")
        atoms_in_res = self.get_atoms(f"resid {resid}")
        old_resname = atoms_in_res[0].resname
        old_atomnames = {a.atomname for a in atoms_in_res}
        new_atomnames = {x[1] for x in types.keys() if x[0] == new_resname}
        old_atomnames_with_matches = set(matching_atoms.keys())
        new_atomnames_with_matches = set(matching_atoms.values())
        # orphans are ones that will be turned into dummies
        orphans = old_atomnames.difference(new_atomnames.union(old_atomnames_with_matches))
        if orphans:
            self.top.parameters.add_dummy_def('DH')
        # newbies are ones that will have to be created from dummies, and bonded to existing ones
        newbies = new_atomnames.difference(old_atomnames.union(new_atomnames_with_matches))
        # first let's add the easy ones (identical + user-matched)
        overlapping = old_atomnames.intersection(new_atomnames)
        masses = gml.guess_element_properties()
        for alch in list(overlapping) + list(old_atomnames_with_matches):
            alch_atom = alch if alch in overlapping else matching_atoms[alch]
            new_type = types[(new_resname, alch_atom)]
            if new_type in self.top.defined_atomtypes:
                new_mass = self.top.parameters.get_subsection('atomtypes').get_entries_by_types(new_type)[0].modifiers[1]
            elif new_type[0] in masses.keys():
                new_mass = masses[new_type[0].upper()][1]
            else:
                new_mass = 1.008
            self.gen_state_b(alch, old_resname, resid, new_type=new_type,
                             new_charge=charges[(new_resname, alch_atom)], new_mass=float(new_mass))
        # then let's turn orphans into dummies
        for alch in orphans:
            self.gen_state_b(alch, old_resname, resid, new_type='DH', new_charge=0.0, new_mass=1.008)
        # then let's add new ones
        for new_alch in newbies:
            last_atom = self.get_atoms(f"resid {resid}")[-1]
            new_type = types[(new_resname, new_alch)]
            if new_type in self.top.defined_atomtypes:
                new_mass = self.top.parameters.get_subsection('atomtypes').get_entries_by_types(new_type)[0].modifiers[
                    1]
            elif new_type[0] in masses.keys():
                new_mass = masses[new_type[0].upper()][1]
            else:
                new_mass = 1.008
            self.add_atom(last_atom.num + 1, "D" + new_alch, atom_type='DH', charge=0.0, resid=resid,
                          resname=old_resname, mass=1.008)
            new_bonds = [bds for bds in bonds[new_resname] if new_alch in bds]
            for new_bond in new_bonds:
                other = [atnm for atnm in new_bond if atnm != new_alch][0]
                if other in newbies:
                    other = 'D' + other
                elif other in matching_atoms.values():
                    other = [k for k, v in matching_atoms.items() if v == other][0]
                try:
                    other_num = self.get_atom(f"resid {resid} and name {other}").num
                except:
                    continue
                else:
                    self.add_bond(other_num, last_atom.num + 1)
            self.gen_state_b("D" + new_alch, old_resname, resid, new_type=new_type,
                             new_charge=charges[(new_resname, new_alch)], new_mass=float(new_mass))

    def gen_state_b(self, atomname: Optional[str] = None, resname: Optional[str] = None,
                    resid: Optional[int] = None, atomtype: Optional[str] = None, new_type: Optional[str] = None,
                    new_charge: Optional[float] = None, new_mass: Optional[float] = None) -> None:
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
                     resid: Optional[int] = None, atomtype: Optional[str] = None) -> None:
        """
        Collapses alchemical A states, making state B
        the new non-alchemical default state A
        :param remove_dummies: bool, whether to remove B-state dummies
        :param atomname: str, name of the selected atom(s) for which state B will be dropped
        :param resname: str, name of the selected residue(s) for which state B will be dropped
        :param resid: int, number of the selected residue(s) for which state B will be dropped
        :param atomtype: str, type of the selected atom(s) for which state B will be dropped
        :return: None
        """
        if not remove_dummies:
            print("WARNING: dropping state A parameters, but keeping dummies (if exist). To remove all atoms with "
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
        self._check_correct()

    def swap_states(self, atomname: Optional[str] = None, resname: Optional[str] = None,
                     resid: Optional[int] = None, atomtype: Optional[str] = None) -> None:
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
            print("WARNING: dropping all state B parameters, but keeping dummies (if exist). To remove all atoms with "
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
        self._check_correct()

    def compare_molecules(self, other: "gml.SectionMol", threshold: float = 0.001) -> None:
        """
        Compares two molecules that should be similar, and prints out differences in their
        respective subsections or parameters
        :param other: gml.SectionMol, the molecule to be compared against
        :param threshold: above this difference in parameters, differences are reported
        :return: None
        """
        self.top.add_ff_params()
        other.top.add_ff_params()
        own_sub = {sub.header for sub in self.subsections if isinstance(sub, gml.SubsectionBonded)}
        other_sub = {sub.header for sub in other.subsections if isinstance(sub, gml.SubsectionBonded)}
        if own_sub.symmetric_difference(other_sub):
            print(f"The two molecules differ by sections: {own_sub.symmetric_difference(other_sub)}")
        if self.natoms != other.natoms:
            raise RuntimeError(f"Can't compare molecules with different numbers of atoms, "
                               f"{self.natoms} and {other.natoms}")
        for sat, oat in zip(self.atoms, other.atoms):
            if abs(sat.charge - oat.charge) > threshold:
                print(f"Charges on atoms {sat} and {oat} differ by more than 0.001")
            if abs(sat.mass - oat.mass) > threshold:
                print(f"Masses on atoms {sat} and {oat} differ by more than 0.001")
            if sat.type != oat.type:
                print(f"Atom types differ for atoms {sat} and {oat}")
        for sub in own_sub.intersection(other_sub):
            owns = self.get_subsection(sub)
            oths = other.get_subsection(sub)
            ownparams = sorted(owns.entries_bonded)
            othparams = sorted(oths.entries_bonded)
            for wnp, thp in zip(ownparams, othparams):
                if wnp.atom_numbers != thp.atom_numbers and wnp.atom_numbers != thp.atom_numbers[::-1]:
                    if (sub == 'dihedrals' and wnp.interaction_type in '42' and thp.interaction_type in '42'
                            and not set(wnp.atom_numbers).symmetric_difference(set(thp.atom_numbers))):
                        pass
                    else:
                        print(f"In section {sub}, topologies diverge between \n{str(wnp).strip()} "
                              f"and \n{str(thp).strip()}")
                elif any([abs(a-b) > threshold for a, b in zip(wnp.params_state_a, thp.params_state_a)]):
                    if sub == 'dihedrals' and wnp.params_state_a[1] == 0 and thp.params_state_a[1] == 0:
                        continue
                    print(f"In section {sub}, parameters have significantly different values "
                          f"between \n{str(wnp).strip()} and \n{str(thp).strip()}")

    def add_atom(self, atom_number: int, atom_name: str, atom_type: str, charge: float = 0.0,
                 resid: Optional[int] = None, resname: Optional[str] = None, mass: Optional[float] = None,
                 print_added: bool = True) -> None:
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
        :param print_added: bool, whether to print the atom being added
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
                q = input("By default, atom will be assigned to residue {}{}. "
                          "Proceed? [y/n]".format(ref_entry.resname, ref_entry.resid))
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
        if print_added:
            print(fstring.format(atom_number, atom_type, resid, resname, atom_name, atom_number, charge, mass).strip())
        new_entry = gml.EntryAtom(fstring.format(atom_number, atom_type, resid, resname, atom_name, atom_number,
                                                 charge, mass), subs_atoms)
        try:
            position = [n for n, a in enumerate(atoms) if isinstance(a, gml.EntryAtom) and a.num == atom_number][0]
        except IndexError:
            try:
                last_atom = [a for a in atoms if isinstance(a, gml.EntryAtom)][-1].num
            except IndexError:
                last_atom = 0
            if atom_number == last_atom + 1:
                atoms.append(new_entry)
            else:
                raise RuntimeError("Last atom number is {}, "
                                   "cannot create atom nr {}".format(last_atom, atom_number))
        else:
            self.offset_numbering(1, atom_number)
            atoms.insert(position, new_entry)
        self.update_dicts()

    def del_atom(self, atom_number: int, del_in_pdb: bool = True, renumber_in_pdb: bool = True,
                 check_after: bool = False) -> None:
        """
        Removes an atom from the topology, as specified using
        topology numbering (1-based)
        :param atom_number: int, atom number in topology
        :param del_in_pdb: bool, whether to also remove in the bound PDB file
        :param renumber_in_pdb: bool, whether to renumber atoms in the bound PDB file
        :param check_after: bool, whether to perform a plain number check after each atom deletion
        :return: None
        """
        if atom_number > self.natoms:
            raise RuntimeError(f"Can't remove atom {atom_number}, molecule only has {self.natoms} atoms")
        if self.top.pdb:
            matched = self._match_pdb_to_top(atom_number) if del_in_pdb else []
        else:
            matched = None
        self._del_atom(atom_number)
        self._del_params(atom_number)
        self.offset_numbering(-1, atom_number)
        self.update_dicts()
        # checking correct numbers, can add more checks in the future
        if check_after:
            self._check_correct()
        if del_in_pdb:
            if self.top.pdb:
                for to_remove in matched:
                    self.top.pdb.delete_atom(to_remove)
                if renumber_in_pdb:
                    self.top.pdb.renumber_atoms()

    def swap_atom(self, atom_number: int, new_position: int, swap_in_pdb: bool = True) -> None:
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

    def _hide_atom(self, old_pos: int, new_pos: int) -> None:
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

    def _return_atom(self, new_pos: int) -> None:
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

    def _del_atom(self, atom_number: int) -> None:
        subsect_atoms = self.get_subsection('atoms')
        chosen = [e for e in subsect_atoms.entries if isinstance(e, gml.EntryAtom) and e.num == atom_number][0]
        subsect_atoms.entries.remove(chosen)

    def _del_params(self, atom_number: int) -> None:
        for subs in [s.header for s in self.subsections if s.header != 'atoms']:
            try:
                subsections = self.get_subsections(subs)
                for subsection in subsections:
                    to_del = []
                    for entry in subsection.entries_bonded:
                        if atom_number in entry.atom_numbers:
                            to_del.append(entry)
                    for entry in to_del:
                        subsection.entries.remove(entry)
            except KeyError:
                pass

    def _check_correct(self) -> None:
        natoms = self.natoms
        for subs in [s.header for s in self.subsections if s.header != 'atoms']:
            try:
                subsection = self.get_subsection(subs)
                for entry in subsection.entries_bonded:
                    if any([e > natoms for e in entry.atom_numbers]):
                        raise RuntimeError(f"Entry {entry} is invalid, only {self.natoms} atoms in the system")
            except KeyError:
                pass

    @property
    def bonds(self) -> list:
        subsection = self.get_subsection('bonds')
        return [entry.atom_numbers for entry in subsection.entries_bonded]

    def add_bond(self, first_atom: int, second_atom: int) -> None:
        """
        This is just an alias for merge_two if bond is intramolecular
        """
        self.merge_two(self, first_atom, second_atom)

    def add_bonds_from_cutoff(self, cutoff: float, cutoff_h: float = None) -> None:
        """
        Goes through all atom pairs in molecule and adds bonds to any pair
        that has a distance lower than the specified cutoff
        :param cutoff: Distance cutoff below which two atoms will be considered bonded (in Angstroms)
        :param cutoff_h: Same but for bonds involving hydrogen atoms
        :return: None
        """
        if not self.top.pdb:
            raise AttributeError("An associated structure is required to run this fn; use Top.add_pdb() to add one")
        if cutoff_h is not None and {a.element for a in self.top.pdb.atoms} == {''}:
            raise AttributeError("No elements are defined to use alternative cutoff for hydrogens. "
                                 "Run Top.pdb.add_elements() to add them automatically")
        if cutoff_h:
            self.top.print("Detected element properties in the structure file. Mind that hydrogens will be identified "
                           "using the last column in the .pdb")
        pdb_atoms_serials = {self._match_pdb_to_top(i)[0] - 1 for i in range(1, self.natoms+1)}
        pdb_atoms = [n for n, a in enumerate(self.top.pdb.atoms) if a.serial in pdb_atoms_serials] # 0-based indices
        pairs = []
        self.top.print("Looking for bonds in the structure...")
        for n, a1 in enumerate(pdb_atoms):
            for a2 in pdb_atoms[n+1:]:
                at1, at2 = self.top.pdb.atoms[a1], self.top.pdb.atoms[a2]
                if cutoff_h is None:
                    if self.top.pdb._atoms_dist_pbc(at1, at2) < cutoff:
                        pairs.append((a1+1, a2+1))
                else:
                    if at1.element.upper() == 'H' or at2.element.upper() == 'H':
                        if self.top.pdb._atoms_dist_pbc(at1, at2) < cutoff_h:
                            pairs.append((a1+1, a2+1))
                    else:
                        if self.top.pdb._atoms_dist_pbc(at1, at2) < cutoff:
                            pairs.append((a1+1, a2+1))
        self.top.print(f"Adding {len(pairs)} bonds to the topology")
        try:
            _ = self.bonds
        except KeyError:
            self.init_subsection('bonds')
        for pair in pairs:
            if pair not in self.bonds and pair[::-1] not in self.bonds:
                self.add_bond(*pair)

    def find_unbound_atoms(self) -> None:
        for a in self.atoms:
            if not self.select_atoms_bonded_to(a.num):
                print(f"Atom {a.num} ({a.atomname}-{a.resname}{a.resid}) has no binding partners")

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

    def merge_two(self, other: "gml.SectionMol", anchor_own: int, anchor_other: int, offset_resid: bool = False) -> None:
        """
        Creates a new bond by either merging two distinct
        molecules (both being part of the same topology)
        or adding a new bond within a single molecule;
        if anchor_own is less than 1, it will simply
        merge two sections without adding any bonds
        :param other: an SectionMol instance, the other molecule that participates in the bond (can be self)
        :param anchor_own: int, number of the atom that will form the new bond in self
        :param anchor_other: int, number of the atom that will form the new bond in other (or self, if other is self)
        :param offset_resid: bool, whether the second molecule's resid should be offset by # of first molecule's residues
        :return: None
        """
        for anch, mol in zip([anchor_own, anchor_other], [self, other]):
            if anch > mol.natoms:
                raise RuntimeError(f"Index {anch} exceeding atom count ({mol.natoms}) in molecule {mol.mol_name}")
        anchor_other = int(anchor_other)
        anchor_own = int(anchor_own)
        if other is not self:
            other.offset_numbering(self.natoms)
            anchor_other += self.natoms
            if offset_resid:
                other.offset_residues(len(self.residues))
            if self.top.pdb:
                self.top.print("WARNING: if merging molecules that are not consecutive, "
                               "make sure to adjust your PDB numbering")
        if anchor_own > 0:
            self._make_bond(anchor_own, anchor_other, other)
        if other is not self:
            self._merge_fields(other)
        if other.mol_name != self.mol_name:
            self.top.sections.remove(other)
            # the stuff below works but is terribly ugly, we need to have API for manipulating content of Top.system
            system_setup = self.top.sections[-1].get_subsection('molecules')
            system_setup.entries = [e for e in system_setup if other.mol_name not in e]
        elif anchor_own < 0:
            self.top.print("WARNING: Merged two molecules with the same names, make sure you adjust "
                           "your [ system ] section manually")

    def merge_molecules(self, other: "gml.SectionMol") -> None:
        """
        Merges two molecules into one section to enable downstream manipulation;
        note there will be fewer molecules defined in the system
        :param other: gml.SectionMol, another molecule definition
        :return: None
        """
        other.offset_numbering(self.natoms)
        self._merge_fields(other)
        self.top.sections.remove(other)
        # the stuff below works but is terribly ugly, we need to have API for manipulating content of Top.system
        system_setup = self.top.sections[-1].get_subsection('molecules')
        system_setup.entries = [e for e in system_setup if other.mol_name not in e]

    def _merge_fields(self, other: "gml.SectionMol") -> None:
        self.top.print('WARNING watch out for #ifdef POSRES keywords that might get misplaced')
        all_sections_own = set(a.header for a in self.subsections if a.header != 'moleculetype')
        all_sections_other = set(a.header for a in other.subsections if a.header != 'moleculetype')
        for subs in list(all_sections_own.union(all_sections_other)):
            # TODO check for a "conditional" attribute of a section
            try:
                subsection_other = other.get_subsection(subs)
                subsection_own = self.get_subsection(subs)
                subsection_own.add_entries([entry for entry in subsection_other if entry])
                subsection_other.section = subsection_own.section
                self.top.print(f"Merging sections {subs} from two molecules")
            except KeyError:
                pass

    def _make_bond(self, atom_own: int, atom_other: int, other: "gml.SectionMol") -> None:
        try:
            _ = other.bonds
        except KeyError:
            other.init_subsection('bonds')
        new_bond = [tuple(sorted([int(atom_own), int(atom_other)]))]
        new_angles = self._generate_angles(other, atom_own, atom_other)
        new_pairs, new_dihedrals = self._generate_14(other, atom_own, atom_other)
        # TODO remove overlapping pairs between new_bond/new_angles and new_pairs for 4- and 5-membered rings
        # or do we really need it?
        if 'pairs' not in [a.header for a in self.subsections] and 'pairs' not in [a.header for a in other.subsections]:
            pairs_needed = False
        else:
            pairs_needed = True
        subs = ['bonds', 'pairs', 'angles', 'dihedrals'] if pairs_needed else ['bonds', 'angles', 'dihedrals']
        news = [new_bond, new_pairs, new_angles, new_dihedrals] if pairs_needed else [new_bond, new_angles, new_dihedrals]
        for sub, entries in zip(subs, news):
            try:
                subsection = self.get_subsection(sub)
            except KeyError:
                self.init_subsection(sub)
                subsection = self.get_subsection(sub)
            subsection.add_entries([gml.EntryBonded(subsection.fstring.format(*entry, subsection.prmtypes[0]), subsection)
                                    for entry in entries])
            try:
                other.get_subsection(sub)
            except KeyError:
                other.init_subsection(sub)

    def remove_bond(self, at1: int, at2: int) -> None:
        """
        Removes a bond between two atoms, and all the associated bonded terms
        :param at1: int, 1-based index of atom 1
        :param at2:  int, 1-based index of atom 2
        :return: None
        """
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

    def _generate_angles(self, other: "gml.SectionMol", atom_own: int, atom_other: int) -> list:
        """
        Generates new angles when an additional bond is formed
        :param other: SectionMol instance, the other molecule that participates in the bond (can be self)
        :param atom_own: int, 1-based index of atom 1
        :param atom_other:  int, 1-based index of atom 2
        :return: list of tuples of indices of atoms that should create the new angles
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
        :param other: SectionMol instance, the other molecule that participates in the bond (can be self)
        :param atom_own: int, 1-based index of atom 1
        :param atom_other: int, 1-based index of atom 4 (forming a 1-4 pair with atom 1)
        :return: two lists containing tuples of indices of atoms that should create (1) the new 1-4 pairs and
        (2) the new dihedrals
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

    def add_ff_params(self, add_section: str = 'all', force_all: bool = False,
                      external_paramsB: Optional["gml.SectionParam"] = None) -> None:
        """
        Looks for FF parameters to be put for every bonded term in the topology,
        then adds them so that they can be explicitly seen/modified
        :param external_paramsB: gml.SectionParam, with this option parameters for the B state can be added
        from another topology
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
                    ssub.add_ff_params(force_all, external_paramsB)

    def find_used_ff_params(self, section: str = 'all') -> list:
        """
        Finds and returns all the FF parameters that are being used by this molecule
        :param section: str, can be a section name, if 'all' then the four standards are used
        :return: list of gml.EntryParam instances
        """
        used_params = []
        if section == 'all':
            subsections_to_add = ['bonds', 'angles', 'dihedrals', 'cmap']  # TODO more options?
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

    def find_missing_ff_params(self, add_section: str = 'all', fix_by_analogy: Optional[dict] = None,
                               fix_B_from_A: bool = False, fix_A_from_B: bool = False, fix_dummy: bool = False,
                               fix_from_external: Optional[str] = None, once: bool = False) -> None:
        """
        Finds FF parameters that are required by the molecule but missing from the parameter definitions
        :param add_section: str, can be a section name, if 'all' then the four standards are used
        :param fix_by_analogy: dict, if defined then missing parameters can get their values reassigned from similar terms
        :param fix_B_from_A: bool, can assign missing parameters for state B directly from state A
        :param fix_A_from_B: bool, can assign missing parameters for state A directly from state B
        :param fix_dummy: bool, if a dummy is missing (atomtype that starts with D), can be automatically added
        :param fix_from_external:
        :param once: bool, only show a missing parameter once for cleaner output
        :return: None
        """
        if add_section == 'all':
            subsections_to_add = ['bonds', 'angles', 'dihedrals']
        else:
            subsections_to_add = [add_section]
        self.printed = []
        missing = self.get_subsection('atoms').check_defined_types()
        if missing and fix_by_analogy:
            for tp in missing:
                if tp in fix_by_analogy.keys():
                    self.top.parameters.clone_type(atomtype=fix_by_analogy[tp], new_type=tp)
                    print(f"Adding atomtype {tp}, copied from atomtype {fix_by_analogy[tp]}")
        elif missing and fix_dummy:
            for tp in missing:
                if tp.startswith('D'):
                    self.top.parameters.add_dummy_def(tp)
                    print(f"Adding dummy atomtype {tp}")
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

    def label_types(self, add_section: str = 'all') -> None:
        """
        Labels all parameters defined in the molecule by atomtype in the comment
        :param add_section: str, if 'all' then applies to all sections, otherwise can be a list of section names
        :return: None
        """
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

    def hydrogen_mass_repartitioning(self, hmass: float = 3.024, methionine=False) -> None:
        """
        Repartitions the masses from heavy atoms to hydrogens to ensure that each
        hydrogen has the desired mass; this enables the use of a 4-fs time step in
        standard MD simulations
        :param hmass: float, desired mass of the hydrogen atom; default is 3.024
        :param methionine: bool, whether to add repartitioning for methionine (might crash in CHARMM otherwise)
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
            if methionine and atom.resname == "MET":
                if atom.atomname == 'SD':
                    atom.mass -= 4.0
                elif atom.atomname == 'CE':
                    atom.mass += 4.0

    def add_posres(self, keyword: Union[str, None] = 'POSRES', value: Union[float, list] = 500.0, selection=None) -> None:
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

    def add_dihres(self, atomlist: list, angle: float, value: float = 500.0, tol: float = 0.0) -> None:
        """
        Adds a dihedral restraint to the molecule's topology
        :param atomlist: list, 1-based indices of the atoms that compose the dihedral
        :param angle: float, value of the dihedral angle
        :param value: float, force constant for the dihedral restraint
        :param tol: float, tolerance
        :return: None
        """
        try:
            sub = self.get_subsection('dihedral_restraints')
        except KeyError:
            content = ['[ dihedral_restraints ]', '; ai   aj   ak   al  type   phi  dphi  kfac']
            sub = gml.SubsectionBonded(content, self)
            self.subsections.append(sub)
        sub.add_entry(gml.EntryBonded(f"{atomlist[0]:6d} {atomlist[1]:6d} {atomlist[2]:6d} {atomlist[3]:6d} 1 {angle:8.3f} {tol:8.3f} {value:8.3f}", sub))

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

    def mutate_protein_residue(self, resid: int, target: str, rtp: Optional[str] = None, mutate_in_pdb: bool = True) -> None:
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
                     ('VAL', 'HG2'): 'HG12', ('VAL', 'HG3'): 'HG13', ('PRO', 'H'): 'HN'}
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
        types, charges, dihedrals, impropers, improper_type, _ = self.parse_rtp(rtp, remember=True)
        # some residue-specific modifications here
        if targ == 'HIS':
            targ = 'HSD' if ('HSD', 'CA') in types.keys() else 'HID'
        elif targ == 'ASH':
            targ = 'ASPP' if ('ASPP', 'CA') in types.keys() else 'ASH'
        elif targ == 'GLH':
            targ = 'GLUP' if ('GLUP', 'CA') in types.keys() else 'GLH'
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
                           'CD': 'CD1', 'HD': 'HD1', 'HD1': ['HD11', 'HE2'], 'HD2': 'HD12', 'HD3': 'HD13', 'H': 'HN'}
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
            self.add_atom(atnum, atom_add, atom_type=types[(targ, atom_add)], charge=charges[(targ, atom_add)],
                          resid=orig.resid, resname=targ, mass=None, print_added=False)
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
        self._check_correct()

    def cleave_protein(self, after_residue: int, rtp: Optional[str] = None, mutate_in_pdb: bool = True) -> None:
        """
        Cuts a protein chain into two, adding typical termini (C- and N-terminal charged ends)
        :param after_residue: int, resid after which the cut will be introduced
        :param rtp: str, a non-standard .rtp file containing the terminal residue definitions
        :param mutate_in_pdb: bool, whether to modify the .pdb file too
        :return: None
        """
        new_c_term = self.get_atom(f'resid {after_residue} and name C')
        new_n_term = self.get_atom(f'resid {after_residue + 1} and name N')
        if mutate_in_pdb:
            if not self.top.pdb:
                print("No .pdb file bound to the topology, use Top.add_pdb() to add one. Only editing topology for now")
            else:
                cterms_pdb = self._match_pdb_to_top(new_c_term.num)
                nterms_pdb = self._match_pdb_to_top(new_n_term.num)
                for ct, nt in zip(cterms_pdb[::-1], nterms_pdb[::-1]):
                    self.top.pdb.make_term('N', nt)
                    self.top.pdb.make_term('C', ct)
        ff = input("Should we use Amber (a) or CHARMM (c) naming for the termini?\n").strip('()').lower()
        while ff not in ['a', 'c']:
            ff = input("Could not understand, please repeat (a) or (c).\n").strip('()').lower()
        if ff == 'a':
            oterm = 'OC'
            htype = 'H'
            otype = 'O2'
        else:
            oterm = 'OT'
            htype = 'HC'
            otype = 'OC'
        self.remove_bond(new_c_term.num, new_n_term.num)
        self.add_atom(new_c_term.num + 1, oterm + '1', otype, 0, new_c_term.resid, new_c_term.resname, 16.0)
        self.add_bond(new_c_term.num, new_c_term.num + 1)
        self.add_atom(new_n_term.num + 2, 'H3', htype, 0, new_n_term.resid, new_n_term.resname, 1.008)
        self.add_bond(new_n_term.num, new_n_term.num + 2)
        self.add_atom(new_n_term.num + 2, 'H2', htype, 0, new_n_term.resid, new_n_term.resname, 1.008)
        self.add_bond(new_n_term.num, new_n_term.num + 2)
        self.atoms[new_c_term.num + 1].atomname = oterm + '2'
        self.atoms[new_c_term.num + 1].type = otype
        self.atoms[new_n_term.num].atomname = 'H1'
        self.atoms[new_n_term.num].type = htype
        rtp = self.find_rtp(rtp) if rtp is None else rtp
        if not 'rtp':
            raise RuntimeError("Failed to locate .rtp files from a local Gromacs installation, please specify"
                               " a custom .rtp file through rtp='/path/to/rtp/file'")
        types, charges, dihedrals, impropers, improper_type, _ = self.parse_rtp(rtp, remember=True)
        if not('N' + new_n_term.resname in types.keys() and 'C' + new_c_term.resname in types.keys()):
            creplaces, cadds, cimpropers = self.parse_tdb(rtp.replace('rtp', 'c.tdb'), 'COO-')
            nterm_res = 'GLY-NH3+' if new_n_term.resname == 'GLY' else 'PRO-NH3+' if new_n_term.resname == 'PRO' else 'NH3+'
            nreplaces, nadds, nimpropers = self.parse_tdb(rtp.replace('rtp', 'n.tdb'), nterm_res)
            for crep in creplaces:
                atom = self.get_atom(f'resid {new_c_term.resid} and name {crep[0]}')
                atom.type = crep[1]
                atom.mass = float(crep[2])
                atom.charge = float(crep[3])
            for nrep in nreplaces:
                atom = self.get_atom(f'resid {new_n_term.resid} and name {nrep[0]}')
                atom.type = nrep[1]
                atom.mass = float(nrep[2])
                atom.charge = float(nrep[3])
            cterm_charge = [ad for ad in cadds if ad[-4].startswith('O')][0][-2]
            self.atoms[new_c_term.num].charge = float(cterm_charge)
            self.atoms[new_c_term.num + 1].charge = float(cterm_charge)
            cterm_type = [ad for ad in cadds if ad[-4].startswith('O')][0][-4]
            self.atoms[new_c_term.num].type = cterm_type
            self.atoms[new_c_term.num + 1].type = cterm_type
            nterm_charge = [ad for ad in nadds if ad[-4].startswith('H')][0][-2]
            self.atoms[new_n_term.num].charge = float(nterm_charge)
            self.atoms[new_n_term.num + 1].charge = float(nterm_charge)
            self.atoms[new_n_term.num + 2].charge = float(nterm_charge)
            nterm_type = [ad for ad in nadds if ad[-4].startswith('H')][0][-4]
            self.atoms[new_n_term.num].type = nterm_type
            self.atoms[new_n_term.num + 1].type = nterm_type
            self.atoms[new_n_term.num + 2].type = nterm_type
            # TODO add impropers
        else:
            ntypes = types['N' + new_n_term.resname]
            ctypes = types['C' + new_c_term.resname]
            ncharges = charges['N' + new_n_term.resname]
            ccharges = charges['C' + new_c_term.resname]
            for catom in self.get_atoms(f'resid {new_c_term.resid}'):
                catom.type = ctypes[catom.atomname]
                catom.charge = ccharges[catom.atomname]
            for natom in self.get_atoms(f'resid {new_n_term.resid}'):
                natom.type = ntypes[natom.atomname]
                natom.charge = ncharges[natom.atomname]
            # TODO add impropers
        self._check_correct()

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
                   self.top.rtp['impropers'], self.top.rtp['bondedtypes'], self.top.rtp['bonds']
        chargedict, typedict = {}, {}
        impropers, dihedrals = {}, {}
        bonds, angles = {}, {}
        bondedtypes = 0
        rtp_cont = [line for line in open(rtp) if not line.strip().startswith(';')]
        resname = None
        reading_atoms = False
        reading_bonds = False
        reading_angles = False
        reading_impropers = False
        reading_dihedrals = False
        reading_bondedtypes = False
        for line in rtp_cont:
            if line.strip().startswith('[') and line.strip().split()[1] not in ['bondedtypes', 'atoms', 'bonds',
                                                                                'exclusions', 'angles', 'dihedrals',
                                                                                'impropers']:
                resname = line.strip().split()[1]
            if line.strip().startswith('[') and line.strip().split()[1] == 'atoms':
                reading_atoms = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'atoms':
                reading_atoms = False
            if line.strip().startswith('[') and line.strip().split()[1] == 'bonds':
                reading_bonds = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'bonds':
                reading_bonds = False
            if line.strip().startswith('[') and line.strip().split()[1] == 'angles':
                reading_angles = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'angles':
                reading_angles = False
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
            if len(line.strip().split()) == 2 and resname is not None and reading_bonds:
                if resname not in bonds.keys():
                    bonds[resname] = []
                bonds[resname].append(line.strip().split())
            # TODO one day we might need angles + bonded entries with parameters?
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
        self.top.rtp['bonds'] = bonds
        return self.top.rtp['typedict'], self.top.rtp['chargedict'], self.top.rtp['dihedrals'], \
            self.top.rtp['impropers'], self.top.rtp['bondedtypes'], self.top.rtp['bonds']

    def parse_tdb(self, tdb: str, terminus: str) -> (list, list, list):
        """
        Reads and parses a terminal database (.tdb) file
        :param tdb: str, an existing .tdb file to read
        :param terminus: str, the terminal residue to read
        :return: a tuple of lists containing replacements, additions, and extra impropers
        """
        replaces, adds, impropers = [], [], []
        tdb_cont = [line for line in open(tdb) if not line.strip().startswith(';')]
        reading = False
        reading_adds = False
        reading_impropers = False
        reading_replaces = False
        for line in tdb_cont:
            if line.strip().startswith('[') and line.strip().split()[1] not in ['replace', 'delete', 'add','impropers']:
                if line.strip().split()[1] == terminus:
                    reading = True
                else:
                    reading = False
            if line.strip().startswith('[') and line.strip().split()[1] == 'add':
                reading_adds = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'add':
                reading_adds = False
            if line.strip().startswith('[') and line.strip().split()[1] == 'impropers':
                reading_impropers = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'impropers':
                reading_impropers = False
            if line.strip().startswith('[') and line.strip().split()[1] == 'replace':
                reading_replaces = True
            if line.strip().startswith('[') and line.strip().split()[1] != 'replace':
                reading_replaces = False
            if len(line.strip().split()) > 3 and reading and reading_replaces:
                replaces.append(line.strip().split())
            if len(line.strip().split()) > 3 and reading and reading_impropers:
                impropers.append(line.strip().split())
            if len(line.strip().split()) > 3 and reading and reading_adds:
                if line.strip().split()[0] in '123456789':
                    adds.append(line.strip().split())
                else:
                    adds[-1].extend(line.strip().split())
        return replaces, adds, impropers

    def update_dicts(self) -> None:
        """
        When atoms are added/removed/transformed, the dictionaries storing atom properties have to be updated
        :return: None
        """
        self.get_subsection('atoms').get_dicts(force_update=True)

    def list_atoms(self):
        """
        Prints all atoms in the molecule
        :return: None
        """
        for atom in self.atoms:
            print(str(atom).strip())

    def list_bonds(self, by_types: bool = False, by_params: bool = False, by_num: bool = False, returning: bool = False):
        """
        Prints or returns (if returning=True) a list of all bonds in the molecule
        :param by_types: bool, if True then atomtypes are used instead of indices
        :param by_params: bool, if True then parameter sets are used instead of indices
        :param by_num: bool, if True then numbers are used instead of indices
        :param returning: bool, if True then the list is returned (otherwise just printed)
        :return: list or None
        """
        return self._list_bonded('bonds', by_types, by_params, by_num, returning)

    def list_angles(self, by_types: bool = False, by_params: bool = False, by_num: bool = False, returning: bool = False):
        return self._list_bonded('angles', by_types, by_params, by_num, returning)

    def list_dihedrals(self, by_types: bool = False, by_params: bool = False, by_num: bool = False, returning: bool = False,
                       interaction_type: Optional[list]=None):
        return self._list_bonded('dihedrals', by_types, by_params, by_num, returning, interaction_type)

    def to_rtp(self, outname: str = 'out.rtp', generate_hdb: Optional[bool] = None, out_hdb: Optional[str] = None,
               set_bonded: bool = True):
        """
        Takes a (single-residue) molecule and makes it into an .rtp entry
        :param outname: str, name for the .rtp file to be produced (or appended to, if exists)
        :param generate_hdb: bool, whether to generate a corresponding .hdb entry (and rename hydrogen atoms)
        :param out_hdb: str, by default outname.replace('rtp', 'hdb') will be used to store the .hdb entry;
        use this option to overwrite this default behavior
        :return: None
        """
        if os.path.exists(outname):
            mode = 'a'
            print('Found an existing .rtp file with requested name, will append to it')
        else:
            mode = 'w'
        if not all([self.atoms[0].resname == a.resname for a in self.atoms]):
            raise RuntimeError(f'Not all atoms in molecule {self.mol_name} have the same residue name, '
                               f'fix this to write an .rtp entry')
        resname = self.atoms[0].resname
        hdbname = outname.replace('rtp', 'hdb') if out_hdb is None else out_hdb
        if generate_hdb is None:
            ans = ''
            while ans.lower() not in 'yntf':
                ans = input(f"Do you also want to generate a hydrogen database (.hdb) entry in file {hdbname}? "
                            f"This will enable automatic generation of hydrogen atoms, but also rename the hydrogens"
                            f"to follow Gromacs' numbering rules. (y/n)")
            generate_hdb = True if ans.lower() in 'yt' else False
        if not generate_hdb:
            replace = {}
        else:
            replace = self.to_hdb(hdbname)
        for a in self.atoms:
            if a.atomname in replace.keys():
                a.atomname = replace[a.atomname]
        atoms = [(a.atomname, a.type, a.charge, n) for n, a in enumerate(self.atoms, 1)]
        # TODO compare against rtp default types
        bonds = self.list_bonds(returning=True, by_params=set_bonded)
        angles = self.list_angles(returning=True, by_params=set_bonded) if set_bonded else []
        dihedrals = self.list_dihedrals(returning=True, interaction_type=["1", "3", "9"], by_params=set_bonded) if set_bonded else []
        impropers = self.list_dihedrals(returning=True, interaction_type=["2", "4"], by_params=set_bonded)
        with open(outname, mode) as outfile:
            outfile.write(f"[ {resname} ]\n")
            outfile.write(f"  [ atoms ]\n")
            for at in atoms:
                outfile.write(f"  {at[0]:6s} {at[1]:6s} {at[2]:8.4f} {at[3]:5d}\n")
            outfile.write(f"  [ bonds ]\n")
            for bd in bonds:
                outfile.write(" ".join(f"{str(i):9s}" for i in bd) + '\n')
            outfile.write(f"  [ angles ]\n" if set_bonded else '')
            for ang in angles:
                outfile.write(" ".join(f"{str(i):9s}" for i in ang) + '\n')
            outfile.write(f"  [ dihedrals ]\n" if set_bonded else '')
            for dih in dihedrals:
                outfile.write(" ".join(f"{str(i):9s}" for i in dih) + '\n')
            if impropers:
                outfile.write(f"  [ impropers ]\n")
                for imp in impropers:
                    outfile.write(" ".join(f"{str(i):9s}" for i in imp) + '\n')

    def to_hdb(self, outname='out.hdb'):
        """
        Takes a (single-residue) molecule and interactively helps you create an .hdb entry
        :param outname: str, name for the .hdb file to be produced (or appended to, if exists)
        :return: dict, hydrogen replacement names
        """
        if os.path.exists(outname):
            mode = 'a'
            print('Found an existing .hdb file with requested name, will append to it')
        else:
            mode = 'w'
        if not all([self.atoms[0].resname == a.resname for a in self.atoms]):
            raise RuntimeError(f'Not all atoms in molecule {self.mol_name} have the same residue name, '
                               f'fix this to write an .hdb entry')
        heavy_with_h = []
        for atom in self.atoms:
            if atom.ish:
                heavy = atom.bound_atoms[0]
                if heavy.num not in [a.num for a in heavy_with_h]:
                    heavy_with_h.append(heavy)
        lines = []
        renamed = {}
        for n, hwh in enumerate(heavy_with_h):
            hs = [a for a in hwh.bound_atoms if a.ish]
            nhs = [a for a in hwh.bound_atoms if not a.ish]
            numh = len(hs)
            nn = chr(65+n) if n <= 25 else chr(64+n//26) + chr(65+n%26)
            hname = f"H{nn}"
            for nn, i in enumerate(range(len(hs)), 1):
                renamed[hs[i].atomname] = hname + str(nn) if len(hs) > 1 else hname
            nhsnh = [a for a in nhs[0].bound_atoms if a.num != hwh.num]
            if numh == 3:
                lines.append((numh, 4, hname, hwh.atomname, nhs[0].atomname, nhsnh[0].atomname))
            elif numh == 2:
                if (hwh.element == "C" and len(nhs) == 1) or hwh.element == "N":
                    lines.append((numh, 3, hname, hwh.atomname, nhs[0].atomname, nhsnh[0].atomname))
                else:
                    lines.append((numh, 6, hname, hwh.atomname, nhs[0].atomname, nhs[1].atomname))
            else:
                assert numh == 1
                if len(nhs) == 3:
                    lines.append((numh, 5, hname, hwh.atomname, nhs[0].atomname, nhs[1].atomname, nhs[2].atomname))
                elif len(nhs) == 2:
                    lines.append((numh, 1, hname, hwh.atomname, nhs[0].atomname, nhs[1].atomname))
                else:
                    lines.append((numh, 2, hname, hwh.atomname, nhs[0].atomname, nhsnh[0].atomname))
        with open(outname, mode) as outfile:
            outfile.write(f"{self.atoms[0].resname}\t{len(heavy_with_h)} \n")
            for line in lines:
                outfile.write("\t".join([str(i) for i in line]) + "\n")
        return renamed

    def write_atomtypes(self, atfile: str = 'atomtypes.dat'):
        """
        Writes (or appends, if exists) atomtypes.dat based on atomtypes defined in the molecule
        :param atfile: str, path to an atomtypes.dat file if has to be appended
        :return: None
        """
        if os.path.exists(atfile):
            mode = 'a'
            print('Found an existing file with requested name, will append to it')
        else:
            mode = 'w'
        element_properties = {"H": (1, 1.008), "O": (8, 15.999), "C": (6, 12.011), "N": (7, 14.007), "S": (16, 32.06),
                              "P": (15, 30.974), "F": (9, 18.998), "X": (0, 0.0)}
        types = set(a.type for a in self.atoms)
        existing = set()
        if mode == 'a':
            existing.update({line.split()[0] for line in open(atfile)})
        to_adds = list(types.difference(existing))
        with open(atfile, mode) as outfile:
            for to_add in to_adds:
                element = to_add[0].upper() if to_add[0].upper() in element_properties.keys() else 'X'
                outfile.write(f"{to_add:8s}      {element_properties[element][1]:12.5f}\n")

    def _list_bonded(self, term: str, by_types: bool = False, by_params: bool = False, by_num: bool = False,
                     returning: bool = False, interaction_type: Optional[list] = None):
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
                    if by_num:
                        appendable = entry.atom_numbers
                    elif not by_types:
                        appendable = entry.atom_names
                    else:
                        appendable = entry.types_state_a
                    returnable.append(appendable + tuple(params))
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

    def copy_atoms_to_vsn(self, selection: Optional[str] = None, type_basename: Optional[str] = None) -> dict:
        """
        Reproduces a Go-Martini-like structure where selected particles get
        a virtual-site "ghost" twin that can have additional interactions
        implemented through nonbond_params
        :param selection: str, the optional selection of atoms that will be duplicated
        :param type_basename: naming convention for the types
        :return: dict matching atom numbers to VSn type
        """
        selected_atoms = self.get_atoms(selection) if selection is not None else self.atoms
        type_basename = f"molecule_{self.top.molecules.index(self)}" if type_basename is None else type_basename
        matchings = {}
        for npart, atom in enumerate(selected_atoms, 1):
            dumtype = f"{type_basename}_{npart}"
            self.add_dummy_def(dumtype)
            self.add_vsn(atom.num, atom.atomname, dumtype)
            matchings[atom.num] = dumtype
        return matchings

    def add_dummy_def(self, dummy_type: str) -> None:
        """
        Adds a dummy type (charge 0, LJ 0) to atomtypes if needed
        :param dummy_type: str, name for the dummy type
        :return: None
        """
        params = self.top.parameters
        params.add_dummy_def(dummy_type)

    def add_vsn(self, atom: int, atomname: str, atomtype: str) -> None:
        """
        Adds a duplication-type virtual site (on top of an atom) at the end of the molecule.
        If the respective section does not exist in the topology, it will be created.
        :param atom: int, serial numbers of the first atom
        :param atomname: str, name of the new atom
        :param atomtype: str, type of the new atom
        :return: None
        """
        try:
            ssect = self.get_subsection('virtual_sitesn')
        except KeyError:
            self.init_subsection('virtual_sitesn')
            ssect = self.get_subsection('virtual_sitesn')
        ref = self.atoms[atom-1]
        self.add_atom(len(self.atoms) + 1, atomname, atomtype, resid=ref.resid, resname=ref.resname, mass=0.0)
        ssect.add_entry(gml.EntryBonded(f'{len(self.atoms)} 1 {atom}\n', ssect))

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
            self.init_subsection('virtual_sites2')
            ssect = self.get_subsection('virtual_sites2')
        ref = self.atoms[atom1-1]
        self.add_atom(len(self.atoms) + 1, atomname, atomtype, resid=ref.resid, resname=ref.resname, mass=0.0)
        ssect.add_entry(gml.EntryBonded(f'{atom1} {atom2} {len(self.atoms)} 1 {fraction}\n', ssect))

    def add_vs3out(self, vsnum: int, atom1: int, atom2: int, atom3: int, a: float, b: float, c: float, atomname: str,
                   atomtype: str):
        """
        Adds an interpolation-type virtual site (between two atoms) just behind the reference (1st) atom.
        If the respective section does not exist in the topology, it will be created.
        :param vsnum: int, serial number of the virtual site
        :param atom1: int, serial numbers of the first atom
        :param atom2: int, serial numbers of the second atom
        :param atom3: int, serial numbers of the second atom
        :param a: float, parameter a of the expression (r𝑖 + 𝑎r_𝑖𝑗 + 𝑏r_𝑖𝑘 + 𝑐(r_𝑖𝑗 × r_𝑖𝑘 ))
        :param b: float, parameter b of the expression (r𝑖 + 𝑎r_𝑖𝑗 + 𝑏r_𝑖𝑘 + 𝑐(r_𝑖𝑗 × r_𝑖𝑘 ))
        :param c: float, parameter c of the expression (r𝑖 + 𝑎r_𝑖𝑗 + 𝑏r_𝑖𝑘 + 𝑐(r_𝑖𝑗 × r_𝑖𝑘 ))
        :param atomname: str, name of the new atom
        :param atomtype: str, type of the new atom
        :return: None
        """
        try:
            ssect = self.get_subsection('virtual_sites3')
        except KeyError:
            self.init_subsection('virtual_sites3')
            ssect = self.get_subsection('virtual_sites3')
        ref = self.atoms[atom1 - 1]
        if any([atom1 > self.natoms, atom2 > self.natoms, atom3 > self.natoms]):
            raise RuntimeError(f"Can't add the vs3, one of the atoms is outside of the molecule {self.mol_name} in "
                               f"the topology; merge your molecules using merge_two")
        self.add_atom(vsnum, atomname, atomtype, resid=ref.resid, resname=ref.resname, mass=0.0)
        ssect.add_entry(gml.EntryBonded(f'{vsnum} {atom1} {atom2} {atom3} 4 {a} {b} {c}\n', ssect))

    def add_constraint(self, atom1: int, atom2: int, distance: float) -> None:
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
            self.init_subsection('constraints')
            ssect = self.get_subsection('constraints')
        ssect.add_entry(gml.EntryBonded(f'{atom1} {atom2} 2 {distance}\n', ssect))

    def make_stateB_dummy(self, resid: int, orig_name: str, dummy_type: str = 'DH') -> None:
        """
        Creates a state B for an atom that is a dummy (will disappear in an alchemical transformation)
        :param resid: int, number of the residue
        :param orig_name: str, name of the atom in that residue that will be modified
        :param dummy_type: str, name of the dummy (DH by default)
        :return: None
        """
        self.add_dummy_def(dummy_type)
        self.gen_state_b(atomname=orig_name, resid=resid, new_type=dummy_type, new_charge=0, new_mass=1.008)

    def solute_tempering(self, temperatures: list, selection: str = None) -> None:
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

    def scale_rest2_bonded(self, gamma: float, selection: str = None) -> None:
        """
        Modifies bonded terms for the REST2 protocol
        :param gamma: float, the scaling parameter from REST2
        :param selection: str or None, can specify the selection to which the scaling should be restricted
        :return: None
        """
        # get a list of atomtypes to clone
        sel = self.get_atoms(selection) if selection is not None else self.atoms
        types = {at.type.strip('y') for at in sel}
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

    def scale_rest2_charges(self, gamma: float, selection: str = None) -> None:
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

    def alchemical_proton(self, resid: int, rtp: Optional[str] = None, b_is_protonated: bool = False) -> None:
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
        if resname not in ['ASP', 'GLU', 'GLUP', 'ASPP', 'LYN', 'LYS']:
            raise RuntimeError("So far only available for residues ASP/ASPP, GLU/GLUP, LYN/LYS")
        mut_dict = {'ASP': 'B', 'GLU': 'J'}
        deprot_dict = {'GLUP': 'GLU', 'ASPP': 'ASP', 'LYS': 'LYN'}
        if resname in ['ASP', 'GLU', 'LYN']:
            self.mutate_protein_residue(resid, mut_dict[resname])
        else:
            for atom in self.get_atoms(f'resid {resid}'):
                atom.resname = deprot_dict[resname]
            resname = deprot_dict[resname]
        atoms = self.get_atoms(f'resid {resid}')
        rtp = self.find_rtp(rtp)
        types, charges, dihedrals, impropers, improper_type, _ = self.parse_rtp(rtp)
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

    def recalc_qtot(self) -> None:
        """
        Puts the qtot (cumulative charge of the molecule) in the comment
        of an atom entry in the topology (done by pdb2gmx by default)
        :return: None
        """
        charge = 0
        for atom in self.atoms:
            charge += atom.charge
            atom.comment = f' ; qtot {charge:.3f}'

    def set_pairs_fudge(self, fudge_LJ=None, fudge_QQ=None, selection=None) -> None:
        """
        Sets explicit 1-4 parameters by taking standard LJ/charges and converting it
        according to [ pairs ] type 2 potential
        :param fudge_LJ: float, scaling factor for epsilon (by default taken from [ defaults ])
        :param fudge_QQ: float, scaling factor for the charge products (by default taken from [ defaults ])
        :param selection: str, will apply the modification to all pairs that are fully within this selection
        :return: None
        """
        if fudge_LJ is None:
            fudge_LJ = self.top.defaults['fudgeLJ']
            self.top.print(f"Setting fudge_LJ to {fudge_LJ} as specified in [ defaults ]")
        if fudge_QQ is None:
            fudge_QQ = self.top.defaults['fudgeQQ']
            self.top.print(f"Setting fudge_QQ to {fudge_QQ} as specified in [ defaults ]")
        if selection is None:
            chosen = []
        else:
            chosen = self.get_atoms(selection)
        pairs_sect = self.get_subsection('pairs')
        all_ats = self.atoms
        for entry in pairs_sect.entries_bonded:
            at1ind, at2ind = entry.atom_numbers
            at1, at2 = all_ats[at1ind-1], all_ats[at2ind-1]
            if selection is not None and (at1 not in chosen or at2 not in chosen):
                continue
            q1, q2 = at1.charge, at2.charge
            s = self.top.parameters.sigma_ij(at1, at2)
            e = self.top.parameters.epsilon_ij(at1, at2)
            entry.interaction_type = '2'
            entry.params_state_a = [fudge_QQ, q1, q2, s, e * fudge_LJ]

    def set_pairs_nb(self, group_a: str, group_b: str, params: dict) -> None:
        """
        Creates the [ pairs_nb ] section between two groups (look up documentation)
        :param group_a: str, selection for group A
        :param group_b: str, selection for group B
        :param params: dict, should contain q_i * q_j, sigma_ij, and epsilon_ij parameters for each consecutive pair
        :return: None
        """
        try:
            subsect = self.get_subsection('pairs_nb')
        except:
            self.subsections.append(gml.SubsectionBonded(['[ pairs_nb ]'], self))
            subsect = self.get_subsection('pairs_nb')
        entries = []
        counter = 0
        for i in [a.num for a in self.get_atoms(group_a)]:
            for j in [a.num for a in self.get_atoms(group_b)]:
                q, s, e = params[counter]
                entries.append(subsect.yield_entry(f" {i} {j} 1 1.0 {q:.12f} {s:.12f} {e:.12f}"))
                counter += 1
        subsect.add_entries(entries)
        print(f"Added {len(entries)} pairs_nb entries")

    def exclude_from_pairs_nb(self) -> None:
        try:
            subsect = self.get_subsection('exclusions')
        except:
            self.subsections.append(gml.SubsectionBonded(['[ exclusions ]'], self))
            subsect = self.get_subsection('exclusions')
        excls = {}
        nb = self.get_subsection('pairs_nb')
        for en in nb.entries_bonded:
            at1, at2 = en.atom_numbers
            if at1 not in excls.keys():
                excls[at1] = []
            if at2 not in excls.keys():
                excls[at2] = []
            excls[at1].append(at2)
            excls[at2].append(at1)
        entries = []
        for atn in sorted(excls.keys()):
            entries.append(subsect.yield_entry(f"{atn} " + " ".join([f"{i}" for i in excls[atn]])))
        subsect.add_entries(entries)


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

    @property
    def atomtypes(self):
        return self.get_subsection('atomtypes')

    @property
    def bondtypes(self):
        return self.get_subsection('bondtypes')

    @property
    def angletypes(self):
        return self.get_subsection('angletypes')

    @property
    def dihedraltypes(self):
        return self.get_subsection('dihedraltypes')

    def _merge(self) -> None:
        """
        If multiple sections (e.g. [ bondtypes ]) are present in the topology,
        this fn merges them into single sections to avoid searching in all instances
        :return: None
        """
        # TODO doesn't always reduce the number of subsections????
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

    def sigma_ij(self, atom_i: "gml.EntryAtom", atom_j: "gml.EntryAtom") -> float:
        """
        A function that will return a correct value of sigma_ij even if custom combination rules are specified
        :param atom_i: gml.EntryAtom, atom i
        :param atom_j: gml.EntryAtom, atom j
        :return: float, the value of sigma in nm
        """
        lor_ber = 0.5 * atom_i.sigma + 0.5 * atom_j.sigma
        try:
            nb = self.get_subsection('nonbond_params')
        except:
            return lor_ber
        else:
            type_i = atom_i.type
            type_j = atom_j.type
            for entry in nb.entries_param:
                if ((entry.types[0] == type_i and entry.types[1] == type_j) or
                        (entry.types[0] == type_j and entry.types[1] == type_i)):
                    return entry.params[0]
            return lor_ber

    def epsilon_ij(self, atom_i: "gml.EntryAtom", atom_j: "gml.EntryAtom") -> float:
        """
        A function that will return a correct value of epsilon_ij even if custom combination rules are specified
        :param atom_i: gml.EntryAtom, atom i
        :param atom_j: gml.EntryAtom, atom j
        :return: float, the value of sigma in kJ/mol
        """
        # TODO only works for Lorentz-Berthelot so far
        lor_ber = (atom_i.epsilon * atom_j.epsilon)**0.5
        try:
            nb = self.get_subsection('nonbond_params')
        except:
            return lor_ber
        else:
            type_i = atom_i.type
            type_j = atom_j.type
            for entry in nb.entries_param:
                if ((entry.types[0] == type_i and entry.types[1] == type_j) or
                        (entry.types[0] == type_j and entry.types[1] == type_i)):
                    return entry.params[1]
            return lor_ber

    def _get_defines(self):
        for sub in self.subsections:
            for entry in [e for e in sub.entries if not isinstance(e, gml.EntryParam)]:
                if entry.content and entry.content[0] == "#define":
                    self.top.defines[entry.content[1]] = entry.content[2:]

    def sort_dihedrals(self) -> None:
        """
        Sorts dihedrals to make sure wildcards are
        moved to the very end of the file
        :return: None
        """
        for sub in self.subsections:
            if 'dihedral' in sub.header:
                sub.sort()  # TODO if two have same periodicity & atoms, remove one with 0-s (PMX)

    def add_dummy_def(self, dummy_type: str) -> None:
        """
        Adds a dummy atom type with a selected name
        :param dummy_type: str, name of the dummy type (e.g. DH or MW)
        :return: None
        """
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

    def fix_zero_periodicity(self) -> None:
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

    def clone_type(self, atomtype: str, prefix: str = 'y', new_type: Optional[str] = None, only_bonded: bool = False) -> None:
        """
        Generates an exact type of a selected atomtype,
        preserving all interactions with other types
        :param atomtype: str, atomtype to be duplicated
        :param prefix: str, new name will be generated as prefix + original atomtype
        :param new_type: str, directly specify the new name (optional)
        :param only_bonded: bool, if the atomtype exists, will only copy bonded terms
        :return: None
        """
        new_atomtype = prefix + atomtype if new_type is None else new_type
        existing_types = self.top.defined_atomtypes
        if new_atomtype in existing_types and not only_bonded:
            raise RuntimeError(f"Type {new_atomtype} already exists, if you want to add bonded terms only, "
                               f"specify only_bonded=True")
        subsect = self.subsections if not only_bonded else [sub for sub in self.subsections
                                                            if not sub.header == "atomtypes"]
        for sub in subsect:
            to_add = []
            for ent in sub:
                if isinstance(ent, gml.EntryParam) and atomtype in ent.types and new_atomtype not in ent.types:
                    to_add.append(ent)
            for entry in to_add:
                newlines = self.gen_clones(entry, atomtype, new_atomtype)
                newentries = [gml.EntryParam(line.replace('\\', ' '), sub, processed=True) for line in newlines]
                sub.add_entries(newentries)
        self.sort_dihedrals()
        self._remove_symm_dupl(new_atomtype)

    def rename_type(self, atomtype: str, prefix: str = 'y', new_type: Optional[str] = None, only_bonded: bool = False) -> None:
        """
        Works like clone_type, but just renames the existing type instead of duplicating it
        :param atomtype: str, atomtype to be duplicated
        :param prefix: str, new name will be generated as prefix + original atomtype
        :param new_type: str, directly specify the new name (optional)
        :param only_bonded: bool, if the atomtype exists, will only copy bonded terms
        :return: None
        """
        new_atomtype = prefix + atomtype if new_type is None else new_type
        subsect = self.subsections if not only_bonded else [sub for sub in self.subsections
                                                            if not sub.header == "atomtypes"]
        for sub in subsect:
            for ent in sub.entries_param:
                for tp in range(len(ent.types)):
                    if ent.types[tp] == atomtype:
                        ent.types = list(ent.types)
                        ent.types[tp] = new_atomtype
                        ent.types = tuple(ent.types)  # that's too silly but typelists should remain tuples <shrugs>

    def clean_unused(self, used_params: list, section: str = 'all') -> None:
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

    def _remove_symm_dupl(self, new_atomtype: str) -> None:
        for sub in self.subsections:
            if 'dihedral' in sub.header:
                sub._remove_symm(new_atomtype)

    def get_opt_dih(self, types: bool = False) -> list:
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and '9' in sub.prmtypes][0]
        return ss.get_opt_dih(types)

    def get_opt_dih_indices(self) -> list:
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and '9' in sub.prmtypes][0]
        return ss.get_opt_dih_indices()

    def set_opt_dih(self, values) -> None:
        ss = [sub for sub in self.subsections if sub.header == 'dihedraltypes' and '9' in sub.prmtypes][0]
        ss.set_opt_dih(values)

    def add_nbfix(self, type1: str, type2: str, mod_sigma: float = 0.0, mod_epsilon: float = 0.0,
                  action_default: str = 'x', new_sigma: Optional[float] = None, new_epsilon: Optional[float] = None,
                  scale_sigma: Optional[float] = None, scale_epsilon: Optional[float] = None) -> None:
        """
        Generates NBFIX entries for the chosen pair of atomtypes by modifying the current
        Lorentz-Berthelot rules, or an existing NBFIX
        :param type1: str, name of the 1st atomtype in the pair
        :param type2: str, name of the 2nd atomtype in the pair
        :param mod_sigma: float, by how much to increase the LJ sigma (in nm)
        :param mod_epsilon: float, by how much to increase the LJ epsilon (in kJ/mol)
        :param new_sigma: float, if specified instead of mod_sigma sets the value directly (in nm)
        :param new_epsilon: float, if specified instead of mod_epsilon sets the value directly (in kJ/mol)
        :param scale_sigma: float, if specified instead of mod_sigma, multiplies instead of adding
        :param scale_epsilon: float, if specified instead of mod_epsilon, multiplies instead of adding
        :param action_default: str, what to do if an NBFIX already exists (check prompt if that happens)
        :return: None
        """
        # TODO has to be faster for Go-models
        tt = time()
        atp = self.get_subsection('atomtypes')
        sigma1, eps1, sigma2, eps2 = [None] * 4
        for entry in atp:
            if isinstance(entry, gml.EntryParam):
                if entry.types[0] == type1:
                    sigma1, eps1 = entry.params
                if entry.types[0] == type2:
                    sigma2, eps2 = entry.params
        # print(time() - tt)
        if sigma1 is None:
            raise KeyError('Type {} was not found in the atomtype definitions'.format(type1))
        if sigma2 is None:
            raise KeyError('Type {} was not found in the atomtype definitions'.format(type2))
        if new_sigma is None:
            if scale_sigma is None:
                new_sig = 0.5 * (sigma1 + sigma2) + mod_sigma
            else:
                new_sig = 0.5 * scale_sigma * (sigma1 + sigma2)
        else:
            new_sig = new_sigma
        if new_epsilon is None:
            if scale_epsilon is None:
                new_eps = (eps1 * eps2) ** 0.5 + mod_epsilon
            else:
                new_eps = scale_epsilon * (eps1 * eps2) ** 0.5
        else:
            new_eps = new_epsilon
        try:
            nbsub = self.get_subsection('nonbond_params')
        except KeyError:
            self.init_subsection('nonbond_params')
            nbsub = self.get_subsection('nonbond_params')
        print('x', new_eps, eps1, eps2)
        comment = ''
        #print(time() - tt)
        # this part is only important if an entry already exists:
        to_remove = []
        for entry in nbsub:
            if isinstance(entry, gml.EntryParam):
                if (entry.types[0], entry.types[1]) in [(type1, type2), (type2, type1)]:
                    action = action_default
                    while action not in 'mrt':
                        action = input("An entry already exists, shall we replace it from scratch (r), modify existing "
                                       "(m) or terminate (t)?")
                    if action == 't':
                        return
                    elif action == 'm':
                        if new_sigma is None:
                            if scale_sigma is None:
                                new_sig = entry.params[0] + mod_sigma
                            else:
                                new_sig = entry.params[0] * scale_sigma
                                mod_sigma = entry.params[0] - new_sig
                        if new_epsilon is None:
                            if scale_epsilon is None:
                                new_eps = entry.params[1] + mod_epsilon
                            else:
                                new_eps = entry.params[1] * scale_epsilon
                                mod_epsilon = entry.params[1] - new_eps
                            print('y', new_eps)
                        comment = entry.comment
                    # if action == 'r' we leave new_sig and new_eps as they are defined above
                    to_remove.append(entry)
        # print(time() - tt)
        for entry in to_remove:
            nbsub.remove_entry(entry)
        entry_line = "{} {} 1 {} {} ; sigma chg by {}, eps chg by {} {}".format(type1, type2, new_sig, new_eps,
                                                                                mod_sigma, mod_epsilon, comment)
        nbsub.add_entry(gml.Subsection.yield_entry(nbsub, entry_line))
        # print(time() - tt)

    def edit_atomtype(self, atomtype: str, mod_sigma: float = 0.0, mod_epsilon: float = 0.0,
                      new_sigma: Optional[float] = None, new_epsilon: Optional[float] = None) -> None:
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
                     atomic_number: Optional[int] = None) -> None:
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
                dupl = [ent for ent in atp.entries if isinstance(ent, gml.EntryParam) and ent.types[0] == atomtype][0]
                atp.entries = [ent for ent in atp.entries if isinstance(ent, gml.EntryParam)
                               and ent.types[0] != atomtype]
                if (abs(float(mass) - float(dupl.modifiers[1])) < 0.000001 and
                        abs(float(sigma) - float(dupl.params[0])) < 0.000001 and
                        abs(float(epsilon) - float(dupl.params[1])) < 0.000001):
                    print(f"Overwriting entry for atomtype {atomtype} with identical parameters (within 1e-6)")
                else:
                    message = f"Overwriting entry for atomtype {atomtype} with "
                    if abs(float(mass) - float(dupl.modifiers[1])) > 0.000001:
                        message += f"mass {mass} (original {dupl.modifiers[1]}), "
                    if abs(float(sigma) - float(dupl.params[0])) > 0.000001:
                        message += f"sigma {sigma} (original {dupl.params[0]}), "
                    if abs(float(epsilon) - float(dupl.params[1])) > 0.000001:
                        message += f"epsilon {epsilon} (original {dupl.params[1]})"
                    print(message.rstrip(", "))
        atp.add_entry(gml.EntryParam(f'{atomtype} {atnum} {mass} 0.0000 A {sigma} {epsilon}', subsection=atp))

    def add_bonded_param(self, types: tuple, params: list, interaction_type: int, action_default: str = 'x') -> None:
        """
        Adds a bonded parameter to a parameter set, creating the subsection if necessary
        :param types: tuple of str, atom types defining the interaction
        :param params: list, parameters of the interaction
        :param interaction_type: int, ID of the interaction type
        :param action_default: str, default action: 'r' for replacing, 't' for skipping, 'a' for duplicating,
        anything else for interactive selection
        :return: None
        """
        subsection_dict = {(2, 1): 'bondtypes', (3, 1): 'angletypes', (4, 1): 'dihedraltypes', (4, 9): 'dihedraltypes',
                           (4, 4): 'dihedraltypes', (4, 2): 'dihedraltypes', (5, 1): 'cmaptypes'}
        try:
            subs = self.get_subsection(subsection_dict[(len(types), interaction_type)])
        except KeyError:
            self.init_subsection(subsection_dict[(len(types), interaction_type)])
            subs = self.subsections[-1]
        matching = [ent for ent in subs.entries_param if ent.types == types or ent.types == types[::-1]]
        if matching:
            action = action_default
            while action not in 'rta':
                action = input("An entry already exists, shall we replace it (r), append a duplicate (a) "
                               "or terminate (t)?")
            if action == 't':
                return
            elif action == 'a':
                pass
            else:
                subs.entries = [ent for ent in subs.entries if ent not in matching]
        if subsection_dict[(len(types), int(interaction_type))] == 'cmaptypes':
            resolution, values = params[0]
            text = f"{' '.join(types)} 1 {resolution} {resolution} {' '.join(values)}"
            subs.add_entry(gml.EntryParam(text, subsection=subs, processed=True, perres=True))
        elif len(params) <= 3:
            subs.add_entry(gml.EntryParam(f'{" ".join(types)} {str(interaction_type)} '
                                          f'{" ".join([str(x) for x in params])}', subsection=subs))
        elif interaction_type == 9:
            assert len(params) % 3 == 0
            for i in range(len(params) // 3):
                param_slice = params[i*3: (i+1)*3]
                subs.add_entry(gml.EntryParam(f'{" ".join(types)} {str(interaction_type)} '
                                              f'{" ".join([str(x) for x in param_slice])}', subsection=subs))

    def gen_clones(self, entry: "gml.EntryParam", atomtype: str, new_atomtype: str) -> list:
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
            lines.append(self.mod_types(entry, mod, new_atomtype, atomtype))
        return lines

    @staticmethod
    def gen_combs(count: int, tuples: int) -> list:
        return list(combinations(range(count), tuples))

    def mod_types(self, entry: "gml.EntryParam", mods: list, new_atomtype: str, atomtype: str) -> str:
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
