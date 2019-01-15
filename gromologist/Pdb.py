from .Entries import EntryAtom


class Pdb:
    def __init__(self, filename, top=None, altloc='A'):
        self.fname = filename
        self.contents = [line.strip() for line in open(self.fname)]
        self.atoms, self.box, self.remarks = self.parse_contents()
        self.top = top
        self.altloc = altloc
        self.atom_format = "ATOM  {:>5d} {:4s}{:1s}{:4s}{:1s}{:>4d}{:1s}   " \
                           "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
        self.cryst_format = "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n"

    def __repr__(self):
        return "PDB file {} with {} atoms".format(self.fname, len(self.atoms))
    
    def add_chains(self, serials=None, chain=None):
        """
        Given a matching Top instance, adds chain identifiers to atoms
        based on the (previously verified) matching between invididual
        molecules defined in Top and consecutive atoms in this Pdb instance.
        Solvent molecules are ommited, and only characters from A to Z
        are supported as valid chain identifiers.
        Optionally, given `serials` and `chain`, one can set all atoms
        with atom numbers in `serials` as belonging to chain `chain`.
        :param serials: iterable of ints, atom numbers to set chain for (default is all)
        :param chain: chain to set for serials (default is use consecutive letters)
        :return: None
        """
        self.check_top()
        if (serials is None and chain is not None) or (serials is not None and chain is None):
            raise ValueError("Both serials and chain have to be specified simultaneously")
        excluded = {'SOL', 'HOH', 'TIP3', 'K', 'NA', 'CL', 'POT'}
        base_char = 65  # 65 is ASCII for "A"
        index = 0
        for mol_name in self.top.system.keys():
            if mol_name.upper() not in excluded:
                n_mols = self.top.system[mol_name]
                mol = self.top.get_molecule(mol_name)
                n_atoms = mol.natoms
                for m in range(n_mols):
                    for a in range(n_atoms):
                        if serials is None and chain is None:
                            self.atoms[index].chain = chr(base_char)
                        else:
                            if self.atoms[index].serial in serials:
                                self.atoms[index].chain = chain
                        index += 1
                    base_char += 1
                    if base_char > 90:
                        return

    def check_top(self, maxwarn=20, fix_pdb=False, fix_top=False):
        if fix_pdb and fix_top:
            raise ValueError("You either want to fix topology or pdb naming")
        if self.top is None:
            raise ValueError("a Top object has not been assigned; molecule info missing")
        index, err = 0, 0
        self.remove_altloc()
        for mol_name in self.top.system.keys():
            mol = self.top.get_molecule(mol_name)
            n_mols = self.top.system[mol_name]
            atom_subsection = mol.get_subsection('atoms')
            atom_entries = [e for e in atom_subsection if isinstance(e, EntryAtom)]
            for m in range(n_mols):
                for n, a in enumerate(atom_entries):
                    rtrn = self.check_mismatch(atom_entries[n], self.atoms[index], mol_name)
                    if rtrn:
                        if fix_pdb:
                            self.atoms[index].atomname = atom_entries[n].atomname
                        elif fix_top:
                            atom_entries[n].atomname = self.atoms[index].atomname
                    err += rtrn
                    index += 1
                    if err > maxwarn:
                        raise RuntimeError("Error: too many warnings")
    
    @staticmethod
    def check_mismatch(atom_entry, atom_instance, mol_name):
        if atom_entry.atomname != atom_instance.atomname:
            print("Atoms {} ({}) in molecule {} topology and {} ({}) in .pdb have "
                  "non-matching names".format(atom_entry.num, atom_entry.atomname, mol_name,
                                              atom_instance.serial, atom_instance.atomname))
            return 1
        return 0
        
    def remove_altloc(self):
        """
        We only keep one of the alternative locations in case
        there is more (by default, A is kept)
        :return: None
        """
        self.atoms = [a for a in self.atoms if a.altloc in [' ', self.altloc]]
    
    def write_atom(self, atom):
        return self.atom_format.format(atom.serial, atom.atomname, atom.altloc, atom.resname, atom.chain, atom.resnum,
                                       atom.insert, atom.x, atom.y, atom.z, atom.occ, atom.beta, atom.element)
    
    def renumber_all(self):
        for n, atom in enumerate(self.atoms):
            atom.serial = n + 1
    
    def match_order_by_top_names(self):
        if self.top is None:
            raise ValueError("a Top object has not been assigned; molecule info missing")
        new_atoms = []
        for mol_name in self.top.system.keys():
            mol = self.top.get_molecule(mol_name)
            n_mols = self.top.system[mol_name]
            atom_subsection = mol.get_subsection('atoms')
            atom_entries = [e for e in atom_subsection if isinstance(e, EntryAtom)]
            for m in range(n_mols):
                for a in atom_entries:
                    pdb_loc = self.select_atoms("resname {} and resid {} and name {}".format(a.resname, a.resid,
                                                                                             a.atomname))
                    if len(pdb_loc) != 1:
                        raise ValueError("Could not proceed; for match-based renumbering, residue numberings "
                                         "have to be consistent between PDB and .top, atom names need to match, "
                                         "and molecules cannot be repeated.\nError encountered when processing"
                                         "residue {} with resid {}, atom name {}".format(a.resname, a.resid,
                                                                                             a.atomname))
                    new_atoms.append(self.atoms[list(pdb_loc)[0]])
        self.atoms = new_atoms
    
    def insert_atom(self, index, base_atom, **kwargs):
        """
        Inserts an atom into the atomlist. The atom is defined by
        providing a base Atom instance and a number of keyworded modifiers,
        e.g. atomname="CA", serial=15, ...
        :param index: int, the new index of the inserted Atom in the Atoms list
        :param base_atom: Atom, instance based on which the new atom will be based
        :param kwargs: Atom attributes to be held by the new atom
        :return: None
        """
        new_atom = Atom(self.write_atom(base_atom))
        for kw in kwargs.keys():
            if kw not in {"serial", "atomname", "resname", "chain", "resnum", "x", "y", "z", "occ", "beta", "element"}:
                raise ValueError("{} is not a valid Atom attribute")
            new_atom.__setattr__(kw, kwargs[kw])
        self.atoms.insert(index, new_atom)

    def select_atoms(self, selection_string):
        protein_selection = "resname ALA ACE CYS ASP ASPP GLU GLUP PHE GLY HIS HID HIE HSD HSE ILE LYS LEU MET " \
                            "NME NMA ASN PRO GLN ARG SER THR VAL TRP"
        dna_selection = "resname DA DG DC DT DA5 DG5 DC5 DT5 DA3 DG3 DC3 DT3"
        rna_selection = "resname RA RG RC RT RA5 RG5 RC5 RT5 RA3 RG3 RC3 RT3"
        solvent_selection = "resname HOH TIP3 SOL K CL NA"
        selection_string = selection_string.replace('solvent', solvent_selection)
        selection_string = selection_string.replace('water', 'resname HOH TIP3 SOL')
        selection_string = selection_string.replace('protein', protein_selection)
        selection_string = selection_string.replace('nucleic', 'dna or rna')
        selection_string = selection_string.replace('dna', dna_selection)
        selection_string = selection_string.replace('rna', rna_selection)
        return sorted(list(self.select_set_atoms(selection_string)))

    def select_set_atoms(self, selection_string):
        """
        
        :param selection_string: str
        :return:
        """
        assert isinstance(selection_string, str)
        parenth_ranges, operators = self.parse_sel_string(selection_string)
        if not parenth_ranges and not operators:
            if selection_string.strip().startswith("not "):
                return self.find_atoms(selection_string.lstrip()[4:], rev=True)
            else:
                return self.find_atoms(selection_string)
        elif parenth_ranges and not operators:
            if selection_string.strip().startswith("not "):
                set_all = {n for n in range(len(self.atoms))}
                return set_all.difference(self.select_set_atoms(selection_string.strip()[4:].strip('()')))
            else:
                return self.select_set_atoms(selection_string.strip().strip('()'))
        else:
            first_op = selection_string[operators[0][0]:operators[0][1]].strip()
            if first_op == "and":
                return self.select_set_atoms(selection_string[:operators[0][0]])\
                    .intersection(self.select_set_atoms(selection_string[operators[0][1]:]))
            elif first_op == "or":
                return self.select_set_atoms(selection_string[:operators[0][0]]) \
                    .union(self.select_set_atoms(selection_string[operators[0][1]:]))
    
    @staticmethod
    def parse_sel_string(selection_string):
        parenth_ranges = []
        operators = []
        opened_parenth = 0
        beginning = 0
        for nc, char in enumerate(selection_string):
            if char == '(':
                opened_parenth += 1
                if beginning == 0:
                    beginning = nc
            elif char == ')':
                opened_parenth -= 1
                end = nc
                if opened_parenth == 0:
                    parenth_ranges.append((beginning, end))
                    beginning = 0
            if opened_parenth < 0:
                raise ValueError("Improper use of parentheses in selection string {}".format(selection_string))
            if selection_string[nc:nc + 5] == " and " and opened_parenth == 0:
                operators.append((nc, nc + 5))
            if selection_string[nc:nc + 4] == " or " and opened_parenth == 0:
                operators.append((nc, nc + 4))
        if opened_parenth != 0:
            raise ValueError("Improper use of parentheses in selection string {}".format(selection_string))
        return parenth_ranges, operators
    
    def find_atoms(self, sel_string, rev=False):
        chosen = []
        keyw = sel_string.split()[0]
        matchings = {"name": "atomname", "resid": "resnum", "resnum": "resnum", "element": "element",
                     "chain": "chain", "resname": "resname", "serial": "serial"}
        try:
            vals = {int(x) for x in sel_string.split()[1:]}
        except ValueError:
            if " to " in sel_string:
                beg = int(sel_string.split()[1])
                end = int(sel_string.split()[3])
                vals = set(range(beg, end+1))
            else:
                vals = set(sel_string.split()[1:])
        for n, a in enumerate(self.atoms):
            if not rev:
                if a.__getattribute__(matchings[keyw]) in vals:
                    chosen.append(n)
            else:
                if a.__getattribute__(matchings[keyw]) not in vals:
                    chosen.append(n)
        return set(chosen)
        
    def parse_contents(self):
        atoms, remarks = [], []
        box = [7.5, 7.5, 7.5, 90, 90, 90]  # generic default, will be overwritten if present
        for line in self.contents:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atoms.append(Atom(line))
            elif line.startswith("CRYST1"):
                box = [float(line[6+9*a:6+9*(a+1)]) for a in range(3)] + \
                      [float(line[33+7*a:33+7*(a+1)]) for a in range(3)]
            elif not line.startswith('TER') and not line.startswith('END'):
                remarks.append(line)
        return atoms, tuple(box), remarks
        
    def fill_beta(self, values, serials=None):
        """
        Enables user to write arbitrary values to the beta field
        of the PDB entry
        :param values: iterable; values to fill
        :param serials: iterable, optional; can be used to specify a subset of atoms
        :return: None
        """
        if any([v > 999 for v in values]):
            print("At least one value is too large to fit into the `beta` field, consider division "
                  "to make values smaller")
        if not serials:
            serials = [a.serial for a in self.atoms]
        if not all([i == j for i, j in zip(serials, sorted(serials))]):
            raise ValueError("Atoms' serial numbers are not sorted, consider running .renumber_all() first")
        if len(values) != len(serials):
            raise ValueError("Lists 'value' and 'serials' have inconsistent sizes: {} and {}".format(len(values),
                                                                                                     len(serials)))
        index = 0
        for atom in self.atoms:
            if atom.serial in serials:
                atom.beta = values[index]
                index += 1

    def save_pdb(self, outname='out.pdb'):
        with open(outname, 'w') as outfile:
            outfile.write(self.cryst_format.format(*self.box))
            for line in self.remarks:
                outfile.write(line.strip() + '\n')
            for atom in self.atoms:
                outfile.write(self.write_atom(atom))
            outfile.write('ENDMDL\n')


class Atom:
    def __init__(self, line):
        self.label = line[:6].strip()
        self.serial = int(line[6:11].strip())
        self.atomname = line[12:16].strip()
        self.altloc = line[16:17]
        self.resname = line[17:21].strip()
        self.chain = line[21:22]
        self.resnum = int(line[22:26].strip())
        self.insert = line[26:27]
        self.x, self.y, self.z = [float(line[30+8*a:30+8*(a+1)]) for a in range(3)]
        self.occ = float(line[54:60].strip())
        self.beta = float(line[60:66].strip())
        try:
            self.element = line[76:78].strip()
        except IndexError:
            name = self.atomname.strip('1234567890')
            if name in 'CHONSP':
                self.element = name[:1]
            else:
                self.element = name[:2]
        
    def __repr__(self):
        chain = self.chain if self.chain != " " else "unspecified"
        return "Atom {} in residue {}{} of chain {}".format(self.atomname, self.resname, self.resnum, chain)
