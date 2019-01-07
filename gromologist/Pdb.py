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
                        self.atoms[index].chain = chr(base_char)
                        index += 1
                    base_char += 1
                    if base_char > 90:
                        return

    def check_top(self):
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
                    if err > 20:
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
    
    def write_line(self, atom):
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
        box = 0, 0, 0
        for line in self.contents:
            if line.startswith('ATOM'):
                atoms.append(Atom(line))
            elif "CRYST" in line:
                box = []
            elif not line.startswith('TER') and not line.startswith('END'):
                remarks.append(line)
        return atoms, box, remarks
        
    def fill_beta(self, values, serials=None):
        pass

    def save_pdb(self):
        path = self.fname.split('/')
        outname = '/'.join(path)
        with open(outname, 'w') as outfile:
            for line in self.contents:
                outfile.write(line)
            outfile.write('\n')


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
