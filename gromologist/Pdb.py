
class Pdb:
    def __init__(self, filename, top=None):
        self.fname = filename
        self.contents = [line.strip() for line in open(self.fname)]
        self.atoms, self.box, self.remarks = self.parse_contents()
        self.top = top
        self.atom_format = "ATOM  {:>5d} {:4s}{:1s}{:4s}{:1s}{:>4d}{:1s}   " \
                           "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
    
    def add_chains(self):
        pass  # TODO
    
    def write_line(self, atom):
        return self.atom_format.format(atom.serial, atom.atomname, atom.altloc, atom.resname, atom.chain, atom.resnum,
                                       atom.insert, atom.x, atom.y, atom.z, atom.occ, atom.beta, atom.element)

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
