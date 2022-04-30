import gromologist as gml
from collections import defaultdict
import math


class Pdb:
    prot_map = {'ALA': 'A', 'CYS': 'C', 'CYX': 'C', 'CYM': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
                'HIS': 'H', 'HIE': 'H', 'HID': 'H', 'HSD': 'H', 'HSE': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
                'TRP': 'W', 'TYR': 'Y', "GLUP": "E", "ASPP": "D"}

    nucl_map = {'DA': "A", 'DG': "G", 'DC': "C", 'DT': "T", 'DA5': "A", 'DG5': "G", 'DC5': "C", 'DT5': "T",
                'DA3': "A", 'DG3': "G", 'DC3': "C", 'DT3': "T", 'RA': "A", 'RG': "G", 'RC': "C", 'RU': "U",
                'RA5': "A", 'RG5': "G", 'RC5': "C", 'RU5': "U", 'RA3': "A", 'RG3': "G", 'RC3': "C", 'RU3': "U",
                'A': "A", 'G': "G", 'C': "C", 'U': "U", 'A5': "A", 'G5': "G", 'C5': "C", 'U5': "U",
                'A3': "A", 'G3': "G", 'C3': "C", 'U3': "U"}

    def __init__(self, filename=None, top=None, altloc='A', qt=False, **kwargs):
        """
        Initializes a Pdb instnace
        :param filename: str, name of the file (.gro or .pdb)
        :param top: str or Top, a topology matching the structure
        :param altloc: str ('A' or 'B'), which altloc to keep if relevant
        :param qt: bool, whether to try to read the PDBQT format (with charges and types)
        :param kwargs: will be passed to Top constructor if top is a string/path
        """
        self.fname = filename
        if self.fname:
            if self.fname.endswith('gro'):
                self.atoms, self.box, self.remarks = self._parse_contents_gro([line.rstrip()
                                                                               for line in open(self.fname)])
            else:
                self.atoms, self.box, self.remarks = self._parse_contents([line.strip() for line in open(self.fname)],
                                                                          qt)
        else:
            self.atoms, self.box, self.remarks = [], 3 * [100] + 3 * [90], []
        self.top = top if not isinstance(top, str) else gml.Top(top, **kwargs)
        if self.top and not self.top.pdb:
            self.top.pdb = self
        if not self.atoms and self.top:
            self.atoms = [Atom.from_top_entry(entry) for mol in self.top.molecules for entry in mol.atoms]
        self.gbox = self._calc_gro_box()
        self.altloc = altloc
        self.conect = {}
        self._atom_format = "ATOM  {:>5d} {:4s}{:1s}{:4s}{:1s}{:>4d}{:1s}   " \
                            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
        self._atom_format_gro = "{:>5d}{:5s}{:>5s}{:>5d}{:8.3f}{:8.3f}{:8.3f}\n"
        self._cryst_format = "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n"

    def __repr__(self):
        return "PDB file {} with {} atoms".format(self.fname, len(self.atoms))

    def add_top(self, top, **kwargs):
        """
        Adds a Top object to the current Pdb object, enabling some
        operations that couple the two
        :param top: str, path to the .top file
        :param kwargs: dict, extra parameters to pass to Top constructor
        :return: None
        """
        self.top = gml.Top(top, **kwargs)
        if self.top and not self.top.pdb:
            self.top.pdb = self

    @classmethod
    def from_selection(cls, pdb, selection):
        """
        Creates a new Pdb instance as a subset of an existing one,
        given a selection that defines a subset of atoms
        :param pdb: Pdb instance, source structure
        :param selection: str, a syntactically correct selection
        :return: Pdb, a subset of the original structure
        """
        # TODO what to do in case of a bound Top object? implement Top.slice()?
        selected_indices = pdb.get_atom_indices(selection)
        new_pdb = Pdb()
        new_pdb.atoms = [atom for n, atom in enumerate(pdb.atoms) if n in selected_indices]
        new_pdb.box = pdb.box
        new_pdb.remarks = pdb.remarks
        new_pdb.altloc = pdb.altloc
        return new_pdb

    @classmethod
    def from_text(cls, text, orig_pdb):
        """
        Reads coordinates from text, useful for reading trajectories
        :param text: str, contents of the pdb/gro
        :param orig_pdb: Pdb instance, contains general information about the structure (box size, remarks etc)
        :return: a new Pdb instance
        """
        new_pdb = Pdb()
        try:
            ftype = orig_pdb.fname[-3:]
        except AttributeError:
            ftype = orig_pdb.sfname[-3:]
        if ftype == 'pdb':
            new_pdb.atoms = cls._parse_contents([line.strip() for line in text.split('\n')], False)
        elif ftype == 'bqt':
            new_pdb.atoms = cls._parse_contents([line.strip() for line in text.split('\n')], True)
        elif ftype == 'gro':
            new_pdb.atoms = cls._parse_contents_gro([line.strip() for line in text.split('\n')])
        new_pdb.box = orig_pdb.box
        new_pdb.remarks = orig_pdb.remarks
        new_pdb.altloc = orig_pdb.altloc
        return new_pdb

    def print_protein_sequence(self, gaps=False):
        """
        Prints protein sequence chain by chain, recognizing amino acids
        as residues that contain a CA atom; unrecognized residues are written as X
        :param gaps: bool, whether to pad gaps (non-continuous numbering) with dashes '-'
        :return: str, the protein sequence
        """
        chains = list({a.chain.strip() for a in self.atoms if a.atomname == 'CA'})
        sequences = []
        if not chains:
            print("No chains in the molecule, run Pdb.add_chains() first")
            return
        for ch in sorted(chains):
            cas = set(self.get_atom_indices(f'name CA and chain {ch}'))
            atoms = [a for n, a in enumerate(self.atoms) if n in cas and a.altloc in [' ', self.altloc]]
            if not gaps:
                sequences.append(''.join([Pdb.prot_map[i.resname] if i.resname in Pdb.prot_map.keys() else 'X'
                                          for i in atoms]))
            else:
                res = [a.resnum for n, a in enumerate(self.atoms) if n in cas and a.altloc in [' ', self.altloc]]
                seq = []
                counter = 0
                for r in range(res[0], res[-1] + 1):
                    if r in res:
                        rname = atoms[counter].resname
                        counter += 1
                        seq.append(Pdb.prot_map[rname] if rname in Pdb.prot_map.keys() else 'X')
                    else:
                        seq.append('-')
                sequences.append(''.join(seq))
        return sequences

    def print_nucleic_sequence(self):
        """
        Prints nucleic acid sequences chain by chain, recognizing nucleotides
        as residues that contain an O4' atom; unrecognized residues are written as X
        :return: str, the nucleic acid sequence
        """
        mapping = defaultdict(lambda: 'X')
        sequences = []
        mapping.update(Pdb.nucl_map)
        chains = list({a.chain.strip() for a in self.atoms if a.atomname == "O4'"})
        if not chains:
            print("No chains in the molecule, run Pdb.add_chains() first")
            return
        for ch in sorted(chains):
            cas = set(self.get_atom_indices(f"name O4' and chain {ch}"))
            atoms = [a for n, a in enumerate(self.atoms) if n in cas]
            sequences.append(''.join([mapping[i.resname] if i.resname in mapping.keys() else 'X' for i in atoms]))
        return sequences

    def find_missing(self):
        """
        Assuming standard naming conventions, finds atoms that are missing
        in protein structures and prints them (heavy atoms only!)
        :return: None
        """
        pro_bb = ['N', 'O', 'C', 'CA']
        pro_sc = {'A': ['CB'], 'C': ['CB', 'SG'], 'D': ['CB', 'CG', 'OD1', 'OD2'], 'E': ['CB', 'CG', 'CD', 'OE1',
                                                                                         'OE2'],
                  'F': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'G': [], 'H': ['CB', 'CG', 'ND1', 'CE1',
                                                                                      'CD2', 'NE2'],
                  'I': ['CB', 'CG1', 'CG2', 'CD'], 'K': ['CB', 'CG', 'CD', 'CE', 'NZ'], 'L': ['CB', 'CG',
                                                                                              'CD1', 'CD2'],
                  'M': ['CB', 'CG', 'SD', 'CE'], 'N': ['CB', 'CG', 'OD1', 'ND2'], 'P': ['CB', 'CG',
                                                                                        'CD'],
                  'Q': ['CB', 'CG', 'CD', 'OE1', 'NE2'], 'R': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                  'S': ['CB', 'OG'], 'T': ['CB', 'CG2', 'OG1'], 'V': ['CB', 'CG1', 'CG2'], 'W': ['CB', 'CG', 'CD1',
                                                                                                 'CD2', 'NE1', 'CE2',
                                                                                                 'CE3', 'CH2', 'CZ2',
                                                                                                 'CZ3'],
                  'Y': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
                        'OH']}
        alt = {'I': ['CD', 'CD1']}  # pretty temporary and non-extensible, need to rethink
        curr_res = ('X', 0)
        atomlist = []
        for at in self.atoms:
            if (at.resname, at.resnum) != curr_res:
                if 'CA' in atomlist:
                    full = set(pro_bb + pro_sc[Pdb.prot_map[curr_res[0]]])
                    if not full.issubset(set(atomlist)):
                        if Pdb.prot_map[curr_res[0]] in alt.keys():
                            modfull = set([alt[Pdb.prot_map[curr_res[0]]][1] if x == alt[Pdb.prot_map[curr_res[0]]][0]
                                           else x for x in full])
                            if not modfull.issubset(set(atomlist)):
                                print(f"atoms {modfull.difference(set(atomlist))} missing from residue {curr_res}")
                        else:
                            print(f"atoms {full.difference(set(atomlist))} missing from residue {curr_res}")
                # TODO implement for nucleic
                curr_res = (at.resname, at.resnum)
                atomlist = []
            atomlist.append(at.atomname)

    def add_chains(self, selection=None, chain=None, offset=0, maxwarn=100, cutoff=10, protein_only=False, nopbc=False):
        """
        Given a matching Top instance, adds chain identifiers to atoms
        based on the (previously verified) matching between invididual
        molecules defined in Top and consecutive atoms in this Pdb instance.
        Solvent molecules are ommited, and only characters from A to Z
        are supported as valid chain identifiers.
        Optionally, given `serials` and `chain`, one can set all atoms
        with atom numbers in `serials` as belonging to chain `chain`.
        :param selection: str, select only a subset of atoms whose chain ID will be set to the value of 'chain'
        :param chain: str, chain to set for the specified selection
        :param offset: int, start chain ordering from letter other than A
        :param maxwarn: int, max number of warnings before an error shows up
        :param cutoff: float, distance threshold (in A) for chain separation if using geometric criteria
        :param protein_only: bool, whether to only add chains to protein residues
        :param nopbc: bool, whether to ignore PBC distances (assumes molecule is whole)
        :return: None
        """
        base_char = 65 + offset  # 65 is ASCII for "A"
        curr_resid = None
        prev_atom = None
        if (selection is None and chain is not None) or (selection is not None and chain is None):
            raise ValueError("Both serials and chain have to be specified simultaneously")
        if selection is not None and chain is not None:
            for atom in self.get_atoms(selection):
                atom.chain = chain
            return
        if not self.top:
            for atom in self.atoms:
                if not prev_atom:
                    prev_atom = atom
                    curr_resid = atom.resnum
                if atom.resnum != curr_resid:
                    if nopbc:
                        dist = self._atoms_dist(atom, prev_atom)
                    else:
                        # TODO fix
                        dist = self._atoms_dist_pbc(atom, prev_atom)
                    prev_atom = atom
                    curr_resid = atom.resnum
                    if dist > cutoff:
                        if (protein_only and atom.resname in Pdb.prot_map.keys()) or not protein_only:
                            base_char += 1
                if base_char > 90:
                    break
                if (protein_only and atom.resname in Pdb.prot_map.keys()) or not protein_only:
                    atom.chain = chr(base_char)
        else:
            self.check_top(maxwarn=maxwarn)
            excluded = {'SOL', 'HOH', 'TIP3', 'WAT', 'OPC', 'K', 'NA', 'CL', 'POT', 'SOD', 'NA+', 'K+', 'CLA'}
            index = 0
            for mol_name in self.top.system.keys():
                if mol_name.upper() not in excluded:
                    n_mols = self.top.system[mol_name]
                    mol = self.top.get_molecule(mol_name)
                    n_atoms = mol.natoms
                    for m in range(n_mols):
                        for a in range(n_atoms):
                            self.atoms[index].chain = chr(base_char)
                            index += 1
                        if (protein_only and self.atoms[index].resname in Pdb.prot_map.keys()) or not protein_only:
                            base_char += 1
                        if base_char > 90:
                            return

    def check_top(self, maxwarn=20, fix_pdb=False, fix_top=False):
        """
        Checks the structure against the associated topology to identify
        and list potential mismatches in atom names
        :param maxwarn: int, maximum number of warnings to print, default is 20
        :param fix_pdb: bool, whether to set names in Pdb using names from the Top
        :param fix_top: bool, whether to set names in Top using names from the Pdb
        :return: None
        """
        if fix_pdb and fix_top:
            raise ValueError("You either want to fix topology or pdb naming")
        if self.top is None:
            raise ValueError("a Top object has not been assigned; molecule info missing")
        index, err = 0, 0
        self._remove_altloc()
        for atom_top, atom_pdb in zip(self.top.atoms, self.atoms):
            index += 1
            try:
                rtrn = self._check_mismatch(atom_top, atom_pdb, atom_top.molname)
            except IndexError:
                raise RuntimeError("Mismatch encountered: PDB has {} atoms while topology "
                                   "has {}".format(len(self.atoms), len(self.top.atoms)))
            if rtrn:
                if fix_pdb:
                    self.atoms[index].atomname = atom_top.atomname
                elif fix_top:
                    atom_top.atomname = self.atoms[index].atomname
            err += rtrn
            index += 1
            if err > maxwarn > -1:
                raise RuntimeError("Error: too many warnings; use maxwarn=N to allow for up to N exceptions,"
                                   "or -1 to allow for any number of them")
        print("Check passed, all names match")

    @staticmethod
    def _check_mismatch(atom_entry, atom_instance, mol_name):
        """
        Checks a Pdb entry against a Top entry to make sure the atom names match
        :param atom_entry: EntryAtom instance from the Top
        :param atom_instance: Atom instance from the Pdb
        :param mol_name: str, name of the molecule to print in the message
        :return: int, 0 if names match or 1 otherwise
        """
        if atom_entry.atomname != atom_instance.atomname or atom_entry.resname != atom_instance.resname:
            print("Atoms {} ({}/{}) in molecule {} topology and {} ({}/{}) in .pdb have "
                  "non-matching names".format(atom_entry.num, atom_entry.atomname, atom_entry.resname, mol_name,
                                              atom_instance.serial, atom_instance.atomname, atom_instance.resname))
            return 1
        return 0

    def print_mols(self):
        """
        Identifies and lists molecules contained in the structure, chain by chain
        :return: None
        """
        # TODO fix
        chains = list({a.chain.strip() for a in self.atoms})
        if not chains:
            print("No chains in the molecule, run Pdb.add_chains() first")

        def identify(atom):
            first = atom.resname
            if first in Pdb.prot_map.keys() or (first[1:] in Pdb.prot_map.keys() and first[0] in 'NC'):
                return 'Protein'
            elif first in Pdb.nucl_map.keys():
                return 'Nucleic'
            else:
                return atom.resname

        mol_list = []
        mol = identify(self.atoms[0])
        chain = self.atoms[0].chain
        res = self.atoms[0].resnum
        for a in self.atoms:
            if mol in ['Protein', 'Nucleic']:
                if a.chain != chain:
                    mol_list.append([mol, 1])
                    chain = a.chain
                    mol = identify(a)
            else:
                if a.resnum != res:
                    if len(mol) > 1 and mol == mol_list[-1][0]:
                        mol_list[-1][1] += 1
                    else:
                        mol_list.append([mol, 1])
                    mol = identify(a)
                    chain = a.chain
            res = a.resnum
        if len(mol) > 1 and mol == mol_list[-1][0] and mol not in ['Protein', 'Nucleic']:
            mol_list[-1][1] += 1
        else:
            mol_list.append([mol, 1])
        for i in mol_list:
            print(f'{i[0]} {i[1]}')

    def _remove_altloc(self):
        """
        Only keeps one of the alternative locations in case
        there is more (by default, A is kept)
        :return: None
        """
        self.atoms = [a for a in self.atoms if a.altloc in [' ', self.altloc]]

    def _write_atom(self, atom, pdb=True):
        """
        Fn to convert an Atom instance to a line compatible with .pdb or .gro formatting
        :param atom: an Atom instance to be written
        :param pdb: bool, whether to write as .pdb (otherwise .gro)
        :return: str, formatted line
        """
        atom.serial %= 100000
        atom.resnum %= 10000
        if atom.occ > 1000:
            atom.occ %= 1000
        if atom.beta > 1000:
            atom.beta %= 1000
        if pdb:
            return self._atom_format.format(atom.serial, atom.atomname, atom.altloc, atom.resname, atom.chain,
                                            atom.resnum, atom.insert, atom.x, atom.y, atom.z, atom.occ, atom.beta,
                                            atom.element)
        else:
            return self._atom_format_gro.format(atom.resnum, atom.resname, atom.atomname, atom.serial, atom.x / 10,
                                                atom.y / 10, atom.z / 10)

    @staticmethod
    def _write_conect(atom, bonded):
        return 'CONECT' + ('{:>5}' * (len(bonded) + 1)).format(atom, *bonded) + '\n'

    def tip3_to_opc(self, offset=0.147722363):
        """
        Converts a 3-point water model in a structure to a 4-point
        by adding a virtual site to each water molecule
        :param offset: a fractional value used to define the position of the virtual site, for OPC equal to 0.147722363
        :return: None
        """
        names = {"TIP3", "WAT", "SOL", "HOH", "TIP"}
        o_names = {"O", "OW"}
        h_names = {"H", "HW"}
        new_atoms = []
        water_mol = {}
        for n, a in enumerate(self.atoms):
            if a.resname in names:
                if a.atomname in o_names:
                    water_mol["O"] = a
                if a.atomname[:-1] in h_names:
                    try:
                        _ = water_mol["H1"]
                    except KeyError:
                        water_mol["H1"] = a
                    else:
                        water_mol["H2"] = a
                if len(water_mol) == 3:
                    new_atoms.append(a)
                    new_atom = Atom(self._write_atom(water_mol["O"]))
                    new_atom.atomname = "MW"
                    new_atom.set_coords([a + offset * (b + c - 2 * a) for a, b, c
                                         in zip(new_atom.coords, water_mol["H1"].coords, water_mol["H2"].coords)])
                    new_atoms.append(new_atom)
                    water_mol = {}
                else:
                    new_atoms.append(a)
            else:
                new_atoms.append(a)
        self.atoms = new_atoms

    def renumber_atoms(self, offset=1, selection=None):
        """
        Consecutively renumbers all atoms in the structure
        :param offset: int, number of the first atom to be set
        :param selection: str, will only renumber atoms matching the selection
        :return: None
        """
        ats = self.atoms if selection is None else self.get_atoms(selection)
        for n, atom in enumerate(ats, offset):
            atom.serial = n

    def renumber_residues(self, offset=1, selection=None):
        """
        Consecutively renumbers all residues in the structure
        :param offset: int, number of the first residue to be set
        :param selection: str, will only renumber residues matching the selection
        :return: None
        """
        count = offset
        ats = self.atoms if selection is None else self.get_atom_indices(selection)
        for n in ats:
            temp = count
            try:
                if self.atoms[n].resnum != self.atoms[n + 1].resnum or self.atoms[n].chain != self.atoms[n + 1].chain:
                    temp = count + 1
            except IndexError:
                pass
            self.atoms[n].resnum = count
            count = temp

    def reposition_atom_from_hook(self, atomsel, hooksel, bondlength, p1_sel=None, p2_sel=None, vector=None):
        """
        Sets coordinates of an atom based on a "hook" atom (one it's bound to), bond length, and a vector defining
        the direction of the bond
        :param atomsel: str, unique selection for the atom being moved
        :param hooksel: str, unique selection for the hook atom
        :param bondlength: float, length of the bond (in A)
        :param p1_sel: str, selection for the 1st atom defining the vector to position the new atom
        :param p2_sel: str, selection for the 2nd atom defining the vector to position the new atom
        :param vector: iterable, defines the vector direction to position the new atom
        :return: None
        """
        coords = self.get_coords()
        if p1_sel is not None and p2_sel is not None:
            p1_xyz = coords[self.get_atom_index(p1_sel)]
            p2_xyz = coords[self.get_atom_index(p2_sel)]
            vec = [x2 - x1 for x1, x2 in zip(p1_xyz, p2_xyz)]
        elif vector is not None:
            vec = vector
        else:
            raise RuntimeError("In repositioning, please use either p1/p2 selections or specify the vector")
        movable = self.get_atom(atomsel)
        hook_xyz = coords[self.get_atom_index(hooksel)]
        vec_len = sum([x ** 2 for x in vec]) ** 0.5
        scale = bondlength / vec_len
        movable.set_coords([h + scale * v for h, v in zip(hook_xyz, vec)])

    def match_order_by_top_names(self, arange=None):
        """
        Whenever PDB atoms have different ordering than .top ones
        but naming is consistent, we can use the ordering from .top
        to reorder PDB atoms.
        This can be done for the entire system if molecules are unique;
        otherwise, range has to be specified (using the Pythonic convention)
        to avoid ambiguities. In that case, matching molecules will only be
        looked for in the specified range.
        :param arange: tuple, start and end point of the modification (end is excluded)
        :return:
        """
        if self.top is None:
            raise ValueError("a Top object has not been assigned; molecule info missing")
        new_atoms = []
        index = 0
        for mol_name in self.top.system.keys():
            mol = self.top.get_molecule(mol_name)
            n_mols = self.top.system[mol_name]
            atom_subsection = mol.get_subsection('atoms')
            atom_entries = [e for e in atom_subsection if isinstance(e, gml.EntryAtom)]
            for _ in arange(n_mols):
                for a in atom_entries:
                    if not arange or arange[0] <= index < arange[1]:
                        pdb_loc = self.get_atom_indices("resname {} and resid {} and name {}".format(a.resname, a.resid,
                                                                                                     a.atomname))
                        pdb_loc = [loc for loc in pdb_loc if not arange or arange[0] <= index < arange[1]]
                        if len(pdb_loc) != 1:
                            raise ValueError("Could not proceed; for match-based renumbering, residue numberings "
                                             "have to be consistent between PDB and .top, atom names need to match, "
                                             "and molecules cannot be repeated.\nError encountered when processing"
                                             "residue {} with resid {}, atom name {}".format(a.resname, a.resid,
                                                                                             a.atomname))
                        new_atoms.append(self.atoms[list(pdb_loc)[0]])
                    index += 1
        self.atoms = new_atoms

    def insert_atom(self, serial, name="HX", hooksel=None, bondlength=None, p1_sel=None, p2_sel=None,
                    vector=None, renumber=True, **kwargs):
        """
        Inserts an atom into the atomlist. The atom is defined by
        providing a base Atom instance and a number of keyworded modifiers,
        e.g. atomname="CA", serial=15, ...
        :param serial: int, the new serial (1-based) of the inserted Atom in the Atoms list
        :param name: str, name of the new atom
        :param hooksel: str, unique selection for the atom the new atom will be bound to
        :param bondlength: float, length of the bond extending from the "hook" atom (in A)
        :param p1_sel: str, selection for the 1st atom defining the vector to position the new atom
        :param p2_sel: str, selection for the 2nd atom defining the vector to position the new atom
        :param vector: iterable, defines the vector direction to position the new atom
        :param renumber: bool, whether to renumber atoms after inserting
        :param kwargs: Atom attributes to be held by the new atom
        :return: None
        """
        base_atom = self.get_atom(hooksel)
        new_atom = Atom(self._write_atom(base_atom))
        for kw in kwargs.keys():
            if kw not in {"serial", "atomname", "resname", "chain", "resnum", "x", "y", "z", "occ", "beta", "element"}:
                raise ValueError("{} is not a valid Atom attribute")
            new_atom.__setattr__(kw, kwargs[kw])
        self.atoms.insert(serial-1, new_atom)
        if renumber:
            self.renumber_atoms()
        atomsel = f'name {name} and resid {base_atom.resnum} and resname {base_atom.resname}'
        if (atomsel is not None and hooksel is not None and bondlength is not None and
                ((p1_sel is not None and p2_sel is not None) or vector is not None)):
            self.reposition_atom_from_hook(atomsel, hooksel, bondlength, p1_sel, p2_sel, vector)

    def delete_atom(self, serial, renumber=False):
        """
        Removes an atom from the Pdb object
        :param serial: int, serial number of the atom to be removed (usually 1-based)
        :param renumber: bool, whether to renumber atoms after removal
        :return: None
        """
        num = [n for n, a in enumerate(self.atoms) if a.serial == serial]
        if len(num) == 0:
            raise ValueError('No atoms with serial number {}'.format(serial))
        elif len(num) > 1:
            raise ValueError('Multiple atoms with serial number {}; consider renumbering'.format(serial))
        atom = self.atoms.pop(num[0])
        print('Entry {} deleted from PDB'.format(str(atom)))
        if renumber:
            self.renumber_atoms()

    def add_elements(self):
        """
        Guesses the element (last element in the .pdb file) based on the atom name,
        useful for e.g. element-based coloring in VMD
        :return: None
        """
        for a in self.atoms:
            if not a.element:
                a.element = [x for x in a.atomname if not x.isdigit()][0]

    def get_atom_indices(self, selection_string):
        """
        Applies a selection to the structure and returns the 0-based indices of selected atoms
        :param selection_string: str, consistent with the selection language syntax (see README)
        :return: list of int, 0-based (!!!) indices of atoms
        """
        sel = gml.SelectionParser(self)
        return sel(selection_string)

    def get_atom_index(self, selection_string):
        """
        Applies a selection to the structure and returns the 0-based index of the selected atom,
        assuming the selection corresponds to a single unique atom (will throw an error otherwise)
        :param selection_string: str, consistent with the selection language syntax (see README)
        :return: int, 0-based (!!!) index of the selected atom
        """
        sel = gml.SelectionParser(self)
        result = sel(selection_string)
        if len(result) > 1:
            raise RuntimeError("Selection {} returned more than one atom: {}".format(selection_string, result))
        elif len(result) < 1:
            raise RuntimeError("Selection {} returned no atoms".format(selection_string, result))
        return result[0]

    def get_atom(self, selection_string):
        """
        Works as get_atom_index, but returns an atom instead of a list of indices
        :param selection_string: str, consistent with the selection language syntax (see README)
        :return: an Atom instance, the selected atom
        """
        return self.atoms[self.get_atom_index(selection_string)]

    def get_atoms(self, selection_string):
        """
        Works as get_atom_indices, but returns a list of atoms instead of a list of indices
        :param selection_string: str, consistent with the selection language syntax (see README)
        :return: a list of Atom instances, the selected atoms
        """
        return [self.atoms[i] for i in self.get_atom_indices(selection_string)]

    def same_residue_as(self, query_iter):
        """
        Broadens the query to all atoms contained in residues from which atoms were selected
        :param query_iter: iterable of int, indices of the atoms in the query
        :return: set, a broadened list of atom indices
        """
        new_list = []
        for atom in query_iter:
            residue, resid = self.atoms[atom].resname, self.atoms[atom].resnum
            matching = [n for n, a in enumerate(self.atoms) if a.resname == residue and a.resnum == resid]
            new_list.extend(matching)
        return set(new_list)

    def within(self, query_iter, threshold, nopbc=False):
        """
        Returns a set of all atoms contained within the specified radius of a selection
        :param query_iter: iterable of int, indices of the atoms in the query
        :param threshold: float, a distance within which atoms will be included
        :param nopbc: bool, whether to include PBC in the distance calculation
        :return: set, a broadened list of atom indices
        """
        new_list = []
        for n, atom in enumerate(self.atoms):
            if nopbc:
                if any([self._atoms_dist(atom, self.atoms[query]) <= threshold for query in query_iter]):
                    new_list.append(n)
            else:
                if any([self._atoms_dist_pbc(atom, self.atoms[query]) <= threshold for query in query_iter]):
                    new_list.append(n)
        return set(new_list)

    @staticmethod
    def _atoms_dist(at1, at2):
        return ((at2.x - at1.x) ** 2 + (at2.y - at1.y) ** 2 + (at2.z - at1.z) ** 2) ** 0.5

    def _atoms_dist_pbc(self, at1, at2):
        if not self.box[3] == self.box[4] == self.box[5] == 90.0:
            a = [self.gbox[0] * 10, self.gbox[3] * 10, self.gbox[4] * 10]
            b = [self.gbox[5] * 10, self.gbox[1] * 10, self.gbox[6] * 10]
            c = [self.gbox[7] * 10, self.gbox[8] * 10, self.gbox[2] * 10]
            d = self._atoms_vec(at1, at2)
            mindist = []
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        mindist.append(sum([(d[x] + a[x] * i + b[x] * j + c[x] * k) ** 2 for x in range(3)]) ** 0.5)
            return min(mindist)
        else:
            return (min([abs(at2.x - at1.x), self.box[0] - abs(at2.x - at1.x)]) ** 2 +
                    min([abs(at2.y - at1.y), self.box[1] - abs(at2.y - at1.y)]) ** 2 +
                    min([abs(at2.z - at1.z), self.box[2] - abs(at2.z - at1.z)]) ** 2) ** 0.5

    @staticmethod
    def _atoms_vec(at1, at2):
        return at2.x - at1.x, at2.y - at1.y, at2.z - at1.z

    def _atoms_vec_pbc(self, at1, at2):
        a = [self.gbox[0] * 10, self.gbox[3] * 10, self.gbox[4] * 10]
        b = [self.gbox[5] * 10, self.gbox[1] * 10, self.gbox[6] * 10]
        c = [self.gbox[7] * 10, self.gbox[8] * 10, self.gbox[2] * 10]
        d = self._atoms_vec(at1, at2)
        vecs = []
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    vecs.append([(d[x] + a[x] * i + b[x] * j + c[x] * k) for x in range(3)])
        return [min([v[0] for v in vecs], key=lambda x: abs(x)),
                min([v[1] for v in vecs], key=lambda x: abs(x)),
                min([v[2] for v in vecs], key=lambda x: abs(x))]

    @staticmethod
    def _parse_contents(contents, qt):
        """
        A parser to extract data from .pdb files
        and convert them to internal parameters
        :param contents: list of str, contents of the .pdb file
        :param qt: bool, whether the .pdb contains extra charge/type columns
        :return: (list of Atom instances, tuple of floats of len 6, list of str)
        """
        atoms, remarks, conect = [], [], {}
        box = [7.5, 7.5, 7.5, 90, 90, 90]  # generic default, will be overwritten if present
        for line in contents:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atoms.append(Atom(line, qt))
            elif line.startswith("CRYST1"):
                box = [float(line[6 + 9 * a:6 + 9 * (a + 1)]) for a in range(3)] + \
                      [float(line[33 + 7 * a:33 + 7 * (a + 1)]) for a in range(3)]
            elif not line.startswith('TER') and not line.startswith('END'):
                remarks.append(line)
            elif line.startswith('CONECT'):
                conect[int(line[6:11].strip())] = []
                for i in [11, 16, 21, 26]:
                    if line[i:i + 5].strip():
                        conect[int(line[6:11].strip())].append(int(line[i:i + 5].strip()))
            if line.startswith('END') or line.startswith('ENDMDL'):
                break
        return atoms, tuple(box), remarks

    def remove_hydrogens(self):
        """
        Uses the standard naming convention (see Pdb.add_elements()
        to identify and remove hydrogen atoms
        :return: None
        """
        self.add_elements()
        new_list = [a for a in self.atoms if a.element != 'H']
        self.atoms = new_list

    def mutate_protein_residue(self, resid, target, chain=''):
        """
        Mutates a chosen residue to a different one (in a standard rotameric state)
        :param resid: int, number of the residue to be mutated
        :param target: str, single-letter code of the residue to be mutated
        :param chain: str, optional name of the chain
        :return: None
        """
        # TODO check not to modify other residues' names
        self.renumber_atoms()
        chstr = 'chain {} and '.format(chain) if chain else ''
        orig = self.get_atom('{}resid {} and name CA'.format(chstr, resid))
        mutant = gml.ProteinMutant(orig.resname, target)
        mutant.check_chiral(self, orig)
        atoms_add, hooks, geo_refs, bond_lengths, _, afters = mutant.atoms_to_add()
        atoms_remove = mutant.atoms_to_remove()
        for at in atoms_remove:
            equivalents = {'OG': 'OG1', 'HG': 'HG1', 'HG1': 'HG11', 'HG2': 'HG12', 'HG3': 'HG13', 'CG': 'CG1',
                           'CD': 'CD1', 'HD': 'HD1', 'HD1': 'HD11', 'HD2': 'HD12', 'HD3': 'HD13'}
            print("Removing atom {} from resid {} in structure".format(at, resid))
            try:
                atnum = self.get_atom_index('{}resid {} and name {}'.format(chstr, resid, at))
            except RuntimeError:
                atnum = self.get_atom_index('{}resid {} and name {}'.format(chstr, resid, equivalents[at]))
            _ = self.atoms.pop(atnum)
        for atom_add, hook, geo_ref, bond_length, aft in zip(atoms_add, hooks, geo_refs, bond_lengths, afters):
            print("Adding atom {} to resid {} in structure".format(atom_add, resid))
            for n in range(len(geo_ref)):
                if isinstance(geo_ref[n], tuple):
                    for i in geo_ref[n]:
                        try:
                            self.get_atom_index('{}resid {} and name {}'.format(chstr, resid, i))
                        except RuntimeError:
                            continue
                        else:
                            geo_ref[n] = i
                            break
            if isinstance(hook, tuple):
                for hk in hook:
                    try:
                        _ = self.get_atom_index('{}resid {} and name {}'.format(chstr, resid, hk))
                    except RuntimeError:
                        continue
                    else:
                        hook = hk
                        break
            hooksel = '{}resid {} and name {}'.format(chstr, resid, hook)
            atomsel = '{}resid {} and name {}'.format(chstr, resid, atom_add)
            aftnr = None
            # TODO check we're inserting correctly + add unittest
            if isinstance(aft, tuple):
                for n, af in enumerate(aft):
                    try:
                        _ = self.get_atom_index('{}resid {} and name {}'.format(chstr, resid, af))
                    except RuntimeError:
                        continue
                    else:
                        aftnr = self.get_atom_index('{}resid {} and name {}'.format(chstr, resid, af))
                        break
                else:
                    if aftnr is None:
                        raise RuntimeError(f"Didn't find any of the atoms: {aft} in residue {resid}")
            else:
                aftnr = self.get_atom_index('{}resid {} and name {}'.format(chstr, resid, aft))
            if len(geo_ref) == 2:
                p1sel = '{}resid {} and name {}'.format(chstr, resid, geo_ref[0])
                p2sel = '{}resid {} and name {}'.format(chstr, resid, geo_ref[1])
                self.insert_atom(aftnr, name=atom_add, hooksel=hooksel,
                                 bondlength=bond_length,
                                 p1_sel=p1sel, p2_sel=p2sel, atomname=atom_add)
            else:
                vec = self._vector(geo_ref, resid, chain)
                self.insert_atom(aftnr, name=atom_add, hooksel=hooksel,
                                 bondlength=bond_length,
                                 vector=vec, atomname=atom_add)

        for atom in self.get_atoms('{}resid {}'.format(chstr, resid)):
            atom.resname = mutant.target_3l
        self.renumber_atoms()

    def _vector(self, atnames, resid, chain, nopbc=False):
        chstr = 'chain {} and '.format(chain) if chain else ''
        atoms = [self.get_atom('{}resid {} and name {}'.format(chstr, resid, at)) for at in atnames]
        if nopbc:
            vecs = [self._atoms_vec(atoms[0], at2) for at2 in atoms[1:]]
        else:
            vecs = [self._atoms_vec_pbc(atoms[0], at2) for at2 in atoms[1:]]
        nv = len(vecs)
        f = -1 if nv == 3 else 1
        return f * sum(v[0] for v in vecs) / nv, f * sum(v[1] for v in vecs) / nv, f * sum(v[2] for v in vecs) / nv

    def add_vs2(self, resid, name1, name2, vsname='V1', fraction=0.5):
        """
        Adds a virtual site (VS) defined by two atoms, interpolating between
        their coordinates
        :param resid: int, number of the residue that will contain the VS
        :param name1: str, 1st atom to use for the construction of the VS
        :param name2: str, 2nd atom to use for the construction of the VS
        :param vsname: str, name of the new virtual site
        :param fraction: float, where to put the VS (0 = on atom 1, 1 = on atom 2, can be interpolated or extrapolated)
        :return: None
        """
        serial = self.get_atoms(f"resid {resid}")[-1].serial + 2
        a1 = self.get_atom(f"resid {resid} and name {name1}")
        a2 = self.get_atom(f"resid {resid} and name {name2}")
        dist = self._atoms_dist_pbc(a1, a2)
        self.insert_atom(serial, name=vsname, hooksel=f"resid {resid} and name {name1}",
                         bondlength=dist * fraction, p1_sel=f"resid {resid} and name {name1}",
                         p2_sel=f"resid {resid} and name {name2}", atomname=vsname)

    def interatomic_dist(self, resid1=1, resid2=2):
        """
        Calculates all distances between atoms in two selected residues
        :param resid1: int, 1st residue to consider
        :param resid2: int, 2nd residue to consider
        :return: list of float, all calculated interatomic distances
        """
        dists = []
        for atom1 in self.get_atoms(f"resid {resid1}"):
            for atom2 in self.get_atoms(f"resid {resid2}"):
                dists.append(self._atoms_dist_pbc(atom1, atom2))
        return dists

    def check_chiral_aa(self, nopbc=False):
        """
        Checks for correct chirality in amino acids, first in the backbone
        and then for chiral side chains
        :param nopbc: bool, whether to ignore PBC information
        :return: None
        """
        prot_atoms = self.get_atoms('name CA and not resname GLY')
        self.check_chiral(prot_atoms, 'N', 'C', 'HA HA1', nopbc=nopbc)
        ile_atoms = self.get_atoms('name CB and resname ILE')
        self.check_chiral(ile_atoms, 'CA', 'CG2', 'CG1', 'side chain chirality', nopbc=nopbc)
        thr_atoms = self.get_atoms('name CB and resname THR')
        self.check_chiral(thr_atoms, 'CA', 'CG1 CG2', 'OG1 OG2', 'side chain chirality', nopbc=nopbc)

    def check_chiral(self, cent_atoms_list, at1, at2, at3, label='backbone chirality', printing=True, nopbc=False):
        """
        Decides on correct or wrong chirality of selected chiral
        centers by calculating selected dihedrals
        :param cent_atoms_list: list of Atom instances, central atoms around which to check chirality
        :param at1: str, name of the 1st surrounding atom
        :param at2: str, name of the 2nd surrounding atom
        :param at3: str, name of the 3rd surrounding atom
        :param label: str, type of chirality (only to display in the warning)
        :param printing: bool, whether to print warnings (default) or return True if all checks are passed
        :param nopbc: bool, whether to ignore PBC when calculating vectors
        :return: None (if printing), bool (if not printing)
        """
        for at in cent_atoms_list:
            resnum, resname, chain = at.resnum, at.resname, at.chain
            chn = ' and chain ' + chain if chain.strip() else ''
            n = self.get_atom(f'name {at1} and resnum {resnum} and resname {resname} {chn}')
            c = self.get_atom(f'name {at2} and resnum {resnum} and resname {resname} {chn}')
            h = self.get_atom(f'name {at3} and resnum {resnum} and resname {resname} {chn}')
            chi = self._get_chirality(at, n, c, h, nopbc)
            if chi < -2.5 or 0 > chi > -1.5:
                if printing:
                    print(f"Check {label} for residue {resname} num {resnum}, looks a bit off")
                else:
                    return False
            elif chi > 0:
                if printing:
                    print(f"Check {label} for residue {resname} num {resnum}, looks like a D-form")
                else:
                    return False
            if not printing:
                return True

    def _get_chirality(self, *atoms, nopbc=False):
        """
        Calculates a dihedral defined by 4 atoms that
        constitute a chiral center
        :param atoms: list of len 4, atoms defining the chiral center
        :param nopbc: bool, whether to ignore PBC when calculating vectors
        :return: float, the dihedral's value
        """
        if nopbc:
            v1 = self._atoms_vec(atoms[1], atoms[0])
            v2 = self._atoms_vec(atoms[0], atoms[2])
            v3 = self._atoms_vec(atoms[2], atoms[3])
        else:
            v1 = self._atoms_vec_pbc(atoms[1], atoms[0])
            v2 = self._atoms_vec_pbc(atoms[0], atoms[2])
            v3 = self._atoms_vec_pbc(atoms[2], atoms[3])
        n1 = self._normalize(self._cross_product(v1, v2))
        n2 = self._normalize(self._cross_product(v2, v3))
        m1 = self._cross_product(n1, self._normalize(v2))
        x = self._scalar_product(n1, n2)
        y = self._scalar_product(m1, n2)
        return math.atan2(y, x)

    @staticmethod
    def _cross_product(v1, v2):
        """
        Calculates a cross product between two vectors
        :param v1: iterable of floats, len 3
        :param v2: iterable of floats, len 3
        :return: list, vector of length 3
        """
        return v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]

    @staticmethod
    def _scalar_product(v1, v2):
        """
        Calculates a dot product between two vectors
        :param v1: iterable of floats, len 3
        :param v2: iterable of floats, len 3
        :return: list, vector of length 3
        """
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

    @staticmethod
    def _normalize(v):
        """
        Normalizes a vector
        :param v: iterable of floats, len 3
        :return: list, vector of length 3
        """
        norm = (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) ** 0.5
        return [a / norm for a in v]

    @staticmethod
    def _parse_contents_gro(contents):
        """
        A parser to extract data from .gro files
        and convert them to internal parameters
        :param contents: list of str, contents of the .gro file
        :return: (list of Atom instances, tuple of floats of len 6, list of str)
        """
        contents = [x for x in contents if x.strip()]
        atoms, remarks = [], []
        header = contents[0]
        remarks.append("TITLE     {}".format(header))
        natoms = int(contents[1])
        for line in contents[2:2 + natoms]:
            atoms.append(Atom.from_gro(line))
        if len(contents[-1].split()) == 3:
            box = [10 * float(x) for x in contents[-1].split()] + [90., 90., 90.]
        elif len(contents[-1].split()) == 9:
            boxline = [float(x) for x in contents[-1].split()]
            assert boxline[3] == boxline[4] == boxline[6] == 0
            box = [0.0] * 6
            box[0] = boxline[0]
            box[-1] = math.atan(boxline[1] / boxline[5]) if boxline[5] != 0 else 90
            box[1] = boxline[1] / math.sin(box[-1] * math.pi/180)
            box[2] = math.sqrt(boxline[7] ** 2 + boxline[8] ** 2 + boxline[2] ** 2)
            box[-2] = math.acos(boxline[7] / box[2])
            box[-3] = math.acos((boxline[8] * math.sin(box[-1] * math.pi/180)) / box[2]
                                + math.cos(box[-2]) * math.cos(box[-1] * math.pi/180))
            box[0], box[1], box[2] = box[0] * 10, box[1] * 10, box[2] * 10
            box[3], box[4], box[5] = round(box[3] * 180/math.pi, 4), round(box[4] * 180/math.pi, 4), round(box[5], 4)
        else:
            raise RuntimeError('Can\'t read box properties')
        return atoms, tuple(box), remarks

    def add_conect(self, cutoff=1.55):
        """
        Adds CONECT entries to the PDB, using the value of
        cutoff to determine molecule connectivity
        :param cutoff: float, cut-off to determine atom bonding patterns
        :return: None
        """
        try:
            import numpy as np
        except ImportError:
            # TODO implement a numpy-free version
            raise RuntimeError("This functionality requires numpy to be present")
        else:
            coords = np.array(self.get_coords())
            for n, atom in enumerate(self.atoms):
                dists = np.linalg.norm(coords - np.array(atom.coords), axis=1)
                selected = list(np.where(np.logical_and(dists < cutoff, dists > 0.1))[0] + 1)
                self.conect[n + 1] = selected

    def fill_beta(self, values, serials=None, smooth=False, ignore_mem=False):
        """
        Enables user to write arbitrary values to the beta field
        of the PDB entry
        :param values: iterable; values to fill
        :param serials: iterable, optional; can be used to specify a subset of atoms
        :param smooth: if float, defines sigma (in Angstrom) for beta-value smoothing
        :param ignore_mem: bool, allows to ignore memory warnings
        :return: None
        """
        if any([v > 999 for v in values]):
            print("At least one value is too large to fit into the `beta` field, consider division "
                  "to make values smaller")
        if not serials:
            serials = [a.serial for a in self.atoms]
        if not all([i == j for i, j in zip(serials, sorted(serials))]):
            raise ValueError("Atoms' serial numbers are not sorted, consider running .renumber_atoms() first")
        if len(values) != len(serials):
            raise ValueError("Lists 'value' and 'serials' have inconsistent sizes: {} and {}".format(len(values),
                                                                                                     len(serials)))
        index = 0
        if not smooth:
            for atom in self.atoms:
                if atom.serial in serials:
                    atom.beta = values[index]
                    index += 1
        else:
            if len(serials) > 10000 and not ignore_mem:
                raise RuntimeError("Try to restrict the number of atoms (e.g. selecting CA only), or you're risking "
                                   "running out of memory. To proceed anyway, run again with ignore_mem=True")
            try:
                import numpy as np
            except ImportError:
                print("Needs numpy for smoothing, try installing it")
            else:
                atomnums = [n for n, atom in enumerate(self.atoms) if atom.serial in serials]
                coords = np.array(self.get_coords())[atomnums]
                for atom in self.atoms:
                    dists = np.linalg.norm(coords - np.array(atom.coords), axis=1)
                    weights = np.exp(-(dists ** 2 / (2 * smooth)))
                    weights /= np.sum(weights)
                    atom.beta = np.sum(values * weights)

    def interpolate_struct(self, other, num_inter, write=False):
        """
        Generates linearly & equally spaced intermediates between two structures,
        the current (self) and another PDB with the same number of atoms
        :param other: Pdb instance, other endpoint to interpolate between
        :param num_inter: int, number of intermediate structures
        :param write: bool, whether to save the intermediates as separate .pdb files (True) or just return them (False)
        :return: None (if write=True) or list of Pdb instances (if write=False)
        """
        inter = []
        try:
            import numpy as np
            from copy import deepcopy
        except ImportError:
            raise RuntimeError("Needs numpy & deepcopy for interpolating, try installing it")
        else:
            self_atoms = np.array(self.get_coords())
            other_atoms = np.array(other.get_coords())
            for i in range(1, num_inter + 1):
                incr = (other_atoms - self_atoms) / (num_inter + 1)
                pdb = deepcopy(self)
                pdb.set_coords(self_atoms + i * incr)
                inter.append(pdb)
        if write:
            for n, struct in enumerate([self] + inter + [other]):
                struct.save_pdb(f'interpolated_structure_{n}.pdb')
            return
        else:
            return [self] + inter + [other]

    def save_pdb(self, outname='out.pdb'):
        """
        Saves the structure in the PDB format
        :param outname: str, name of the file being produced
        :return: None
        """
        with open(outname, 'w') as outfile:
            outfile.write(self._cryst_format.format(*self.box))
            for line in self.remarks:
                outfile.write(line.strip() + '\n')
            for atom in self.atoms:
                outfile.write(self._write_atom(atom))
            for conect in self.conect.keys():
                outfile.write(self._write_conect(conect, self.conect[conect]))
            outfile.write('ENDMDL\n')

    def save_gro(self, outname='out.gro'):
        """
        Saves the structure in the GRO format
        :param outname: str, name of the file being produced
        :return: None
        """
        # TODO vectors can be misaligned when converting gro to pdb
        with open(outname, 'w') as outfile:
            outfile.write("written by gromologist\n{}\n".format(len(self.atoms)))
            for atom in self.atoms:
                outfile.write(self._write_atom(atom, pdb=False))
            gbox = self._calc_gro_box()
            outfile.write((9 * "{:10.5f}" + "\n").format(*gbox))

    def _calc_gro_box(self):
        if self.box[3] == self.box[4] == self.box[5] == 90.0:
            return "{:10.5f}{:10.5f}{:10.5f}\n".format(*[x / 10 for x in self.box[:3]])
        else:
            gbox = 9 * [0.0]
            conv = math.pi / 180
            gbox[0] = self.box[0] / 10
            gbox[1] = self.box[1] / 10 * math.sin(self.box[5] * conv)
            gbox[7] = self.box[2] / 10 * math.cos(self.box[4] * conv)
            gbox[8] = self.box[2] / 10 * (math.cos(self.box[3] * conv) - math.cos(self.box[4] * conv) * math.cos(
                self.box[5] * conv)) / math.sin(self.box[5] * conv)
            gbox[2] = math.sqrt(self.box[2] * self.box[2] / 100 - gbox[7] * gbox[7] - gbox[8] * gbox[8])
            gbox[5] = self.box[1] / 10 * math.cos(self.box[5] * conv)
            return gbox

    def get_coords(self, subset=None):
        """
        Returns all atomic coordinates
        :param subset: list of int, 0-based indices of atoms for which coordinates should be retrieved
        :return: list of list of float, coordinates of all atoms
        """
        if subset:
            return [[a.x, a.y, a.z] for a in [self.atoms[q] for q in subset]]
        else:
            return [[a.x, a.y, a.z] for a in self.atoms]

    def set_coords(self, new_coords):
        """
        Sets all atomic coordinates
        :param new_coords: list of list of int, new coordinates to be set
        :return: None
        """
        assert len(new_coords) == len(self.atoms)
        for atom, coords in zip(self.atoms, new_coords):
            assert len(coords) == 3
            atom.x, atom.y, atom.z = coords


class Atom:
    def __init__(self, line, qt=False):
        """
        Represents a single atom contained in the structure file
        :param line: str, line from the structure file
        :param qt: bool, whether there are extra charge (q) and type (t) fields at the end
        """
        self.label = line[:6].strip()
        try:
            self.serial = int(line[6:11].strip())
        except ValueError:
            try:
                self.serial = int(line[6:11].strip(), 16)
            except ValueError:
                raise RuntimeError(f"Cannot interpret atom number {line[6:11].strip()} as decimal or hexadecimal")
        self.atomname = line[12:16].strip()
        self.altloc = line[16:17]
        self.resname = line[17:21].strip()
        self.chain = line[21:22]
        try:
            self.resnum = int(line[22:26].strip())
        except ValueError:
            try:
                self.resnum = int(line[22:26].strip(), 16)
            except ValueError:
                raise RuntimeError(f"Cannot interpret residue number {line[22:26].strip()} as decimal or hexadecimal")
        self.insert = line[26:27]
        self.x, self.y, self.z = [float(line[30 + 8 * a:30 + 8 * (a + 1)]) for a in range(3)]
        if not qt:
            self.occ = float(line[54:60].strip())
            self.beta = float(line[60:66].strip())
            self.q = 0.0
            self.type = 'X'
        else:
            self.occ = 1.0
            self.beta = 0.0
            self.q = float(line[54:62].strip())
            self.type = str(line[62:].strip())
        try:
            self.element = line[76:78].strip()
        except IndexError:
            name = self.atomname.strip('1234567890')
            if name in 'CHONSP':
                self.element = name[:1]
            else:
                self.element = name[:2]

    @classmethod
    def from_gro(cls, line):
        """
        Reads fields from a line formatted according
        to the .gro format specification
        :param line: str, line from a .gro file
        :return: an Atom instance
        """
        data = "ATOM  {:>5d} {:4s} {:4s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00"
        resnum = int(line[:5].strip()) % 10000
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnum = int(line[15:20].strip())
        x, y, z = [float(line[20 + 8 * i:20 + 8 * (i + 1)].strip()) * 10 for i in range(3)]
        return cls(data.format(atomnum, atomname[:4], resname[:4], resnum, x, y, z))

    @classmethod
    def from_top_entry(cls, entry):
        """
        Creates a dummy Atom instance (no coordinates)
        based on atom information provided by the .top file
        :param entry: a gromologist.EntryAtom instance
        :return: an Atom instance
        """
        data = "ATOM  {:>5d} {:4s} {:4s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00"
        resnum = entry.resid
        resname = entry.resname
        atomname = entry.atomname
        atomnum = entry.num
        x, y, z = [0, 0, 0]
        return cls(data.format(atomnum, atomname, resname, resnum, x, y, z))

    @property
    def coords(self):
        """
        The coordinates of a given atom
        :return: list of float
        """
        return [self.x, self.y, self.z]

    def set_coords(self, coords):
        """
        A setter for Atom coordinates
        :param coords: iterable of len 3
        :return: None
        """
        self.x, self.y, self.z = coords

    def __repr__(self):
        chain = self.chain if self.chain != " " else "unspecified"
        return "Atom {} in residue {}{} of chain {}".format(self.atomname, self.resname, self.resnum, chain)


class Traj(Pdb):
    """
    Still not implemented correctly!
    """

    def __init__(self, structfile=None, compressed_traj=None, top=None, array=None, altloc='A', **kwargs):
        super().__init__(structfile, top, altloc, **kwargs)
        self.sfname = structfile
        self.ctfname = compressed_traj
        self.array = array
        self.box_traj = []
        if self.sfname is None and self.ctfname is None:
            self.coords = [[[]]]
        elif self.sfname is None and self.sfname.endswith('.pdb'):
            self.coords = self.get_coords_from_pdb()
        elif self.sfname is None and self.sfname.endswith('.gro'):
            self.coords = self.get_coords_from_gro()
        if self.ctfname is not None and self.ctfname.endswith('.xtc') and self.sfname is not None:
            import mdtraj as md
            traj = md.load(self.ctfname, top=self.sfname)
            self.coords = [frame.xyz[0] for frame in traj]
        if self.array and self.coords:
            self.coords = self.array
        self.top = top if not isinstance(top, str) else gml.Top(top, **kwargs)
        if self.top and not self.top.pdb:
            self.top.pdb = self
        self.altloc = altloc
        self.nframes = len(self.coords)
        self._atom_format = "ATOM  {:>5d} {:4s}{:1s}{:4s}{:1s}{:>4d}{:1s}   " \
                            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
        self._cryst_format = "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n"

    def __repr__(self):
        return "trajectory file {} with {} frames and {} atoms".format(self.fname, self.nframes, len(self.atoms))

    def get_coords_from_pdb(self):
        coords = []
        content = [line for line in open(self.sfname)]
        term_lines = [0] + [n for n, line in enumerate(content)
                            if line.startswith('END') or line.startswith('ENDMDL')]
        for i, j in zip(term_lines[:-1], term_lines[1:]):
            frame = ''.join(content[i:j])
            struct = Pdb.from_text(frame, self)
            coords.append(struct.get_coords())
            # TODO extract box size & append to self.box_traj
        return coords

    def get_coords_from_gro(self):
        coords = []
        content = [line for line in open(self.sfname)]
        natoms = int(content[1].strip())
        per_frame = natoms + 3
        term_lines = range(0, len(content) + 1, per_frame)
        for i, j in zip(term_lines[:-1], term_lines[1:]):
            struct = Pdb.from_text('\n'.join(content[i:j]), self)
            coords.append(struct.get_coords())
            # TODO extract box size & append to self.box_traj
        return coords

    def __str__(self):
        if not self.box_traj:
            text = self._cryst_format.format(*self.box)
        else:
            raise NotImplementedError("Variable box size not yet implemented")
        for frame in range(len(self.coords)):
            text = text + 'MODEL     {:>4d}\n'.format(frame + 1)
            self.set_coords(self.coords[frame])
            for atom in self.atoms:
                text = text + self._write_atom(atom)
            text = text + 'ENDMDL\n'
        return text

    def save_traj_as_pdb(self, filename='traj.pdb'):
        with open(filename, 'w') as outfile:
            outfile.write(str(self))
