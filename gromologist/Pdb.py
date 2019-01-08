from .Entries import EntryAtom
import warnings


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

    def check_top(self, maxwarn=20):
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
                    err += self.check_mismatch(atom_entries[n], self.atoms[index], mol_name)
                    index += 1
                    if err > maxwarn:
                        raise RuntimeError("Error: too many warnings")
    
    @staticmethod
    def check_mismatch(atom_entry, atom_instance, mol_name):
        if atom_entry.atomname != atom_instance.atomname:
            warnings.warn("Atoms {} ({}) in molecule {} topology and {} ({}) in .pdb have "
                          "non-matching names".format(atom_entry.num, atom_entry.atomname, mol_name,
                                                      atom_instance.serial, atom_instance.atomname), Warning)
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

    def find_line(self, name, res, resname):
        lines = [x for x in range(len(self.contents)) if self.contents[x].startswith("ATOM")
                 and int(self.contents[x][22:26]) == res and self.contents[x][12:16].strip() == name
                 and self.contents[x][17:20].strip() == resname]
        if len(lines) > 1:
            raise RuntimeError("found more than one line fulfilling criterion:\n{}\n"
                               "Consider renumbering DNA residues in your PDB file to avoid duplicates"
                               .format(" ".join([self.contents[a] for a in lines])))
        elif len(lines) == 0:
            raise RuntimeError("name {} and resnum {} not found".format(name, res))
        return lines[0]
    
    def parse_contents(self):
        atoms, remarks = [], []
        box = [7.5, 7.5, 7.5, 90, 90, 90]  # generic default, will be overwritten if present
        for line in self.contents:
            if line.startswith('ATOM'):  # TODO take care of HETATM
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
        if any([v > 9999 for v in values]):
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
            self.element = line[76:78]
        except IndexError:
            name = self.atomname.strip('1234567890')
            if name in 'CHONSP':
                self.element = name[:1]
            else:
                self.element = name[:2]
        
    def __repr__(self):
        chain = self.chain if self.chain != " " else "unspecified"
        return "Atom {} in residue {}{} of chain {}".format(self.atomname, self.resname, self.resnum, chain)
