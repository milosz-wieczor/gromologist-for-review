import os
from .Section import *


class Top:
    def __init__(self, filename, gmx_dir='/usr/share/gromacs/top/'):
        """
        A class to represent and contain the Gromacs topology file and provide
        tools for editing topology elements
        :param filename: str, path to the .top file
        :param gmx_dir: str, Gromacs FF directory
        """
        self.fname = filename
        self.top = self.fname.split('/')[-1]
        self.dir = '/' + '/'.join(self.fname.split('/')[:-1])
        self.contents = open(self.fname).readlines()
        self.gromacs_dir = gmx_dir
        self.include_all()
        self.sections = []
        self.parse_sections()
        self.system, self.charge, self.natoms = self.read_system_properties
    
    def include_all(self):
        """
        includes all .itp files in the .top file to facilitate processing
        :return: None
        """
        # TODO do we want to add #ifdef support + .mdp processing?
        lines = [i for i in range(len(self.contents)) if self.contents[i].startswith("#include")]
        while len(lines) > 0:
            lnum = lines[0]
            to_include = self.find_in_path(self.contents.pop(lnum).split()[1].strip('"\''))
            contents = open(to_include).readlines()
            self.contents[lnum:lnum] = contents
            lines = [i for i in range(len(self.contents)) if self.contents[i].startswith("#include")]
    
    def find_in_path(self, filename):
        """
        looks for a file to be included in either the current directory
        or in Gromacs directories (as given by user), in order to
        include all .itp files in a single .top file
        :param filename: str, name of the file to be searched for
        :return: None
        """
        pref = ''
        if filename in os.listdir(self.dir):
            return self.dir + '/' + filename
        else:
            if '/' in filename:
                pref = '/'.join(filename.split('/')[:-1])
            suff = filename.split('/')[-1]
            if pref.startswith('/') or pref.startswith('./'):
                first = pref.split('/')[1]
            else:
                first = pref.split('/')[0]
            if first in os.listdir(self.dir) and suff in os.listdir(self.dir + '/' + pref):
                return self.dir + '/' + pref + '/' + suff
            elif suff in os.listdir(self.gromacs_dir + '/' + pref):
                return self.gromacs_dir + '/' + pref + '/' + suff
        raise ValueError('file {} not found in neither local nor Gromacs directory'.format(filename))
    
    def parse_sections(self):
        """
        Cuts the content in sections as defined by the position
        of special headers, and builds the self.sections list
        :return: None
        """
        special_sections = {'defaults', 'moleculetype', 'system'}
        special_lines = [n for n, l in enumerate(self.contents) if l.strip().strip('[]').strip() in special_sections]
        special_lines.append(len(self.contents) - 1)
        for beg, end in zip(special_lines[:-1], special_lines[1:]):
            self.sections.append(self.yield_sec(self.contents[beg:end]))
    
    def yield_sec(self, content):
        """
        Chooses which class (Section or derived classes)
        should be used for the particular set of entries
        :param content: list of str, slice of self.content
        :return: Section (or its subclass) instance
        """
        if 'defaults' in content[0]:
            return SectionParam(content, self)
        elif 'moleculetype' in content[0]:
            return SectionMol(content, self)
        else:
            return Section(content, self)
        
    def read_system_properties(self):
        system_subsection = [s.get_subsection('molecules') for s in self.sections if isinstance(s, Section)
                             and 'molecules' in [ss.header for ss in s.subsections]]
        molecules = {}
        natoms, charge = 0, 0
        if len(system_subsection) == 0:
            raise KeyError
        elif len(system_subsection) > 1:
            print("Section 'molecules' not present in the topology")
            return None, None, None
        for e in system_subsection:
            if not e.startswith(';'):
                molecules[e.split()[0]] = int(e.split()[1])
        for mol in molecules.keys():
            sect_mol = self.get_molecule(mol)
            natoms += molecules[mol] * sect_mol.natoms
            charge += molecules[mol] * sect_mol.charge
        return molecules, charge, natoms
        
    def get_molecule(self, mol_name):
        """
        Finds a molecule (SectionMol instance) whose mol_name
        matches the query name
        :param mol_name: name of the molecule
        :return: SectionMol instance
        """
        mol = [s for s in self.sections if isinstance(s, SectionMol) and s.mol_name == mol_name]
        if len(mol) == 0:
            raise KeyError
        elif len(mol) > 1:
            raise RuntimeError("Molecule {} is duplicated in topology".format(mol_name))
        return mol[0]
    
    def check_pdb(self, pdb_object):
        """
        Compares the topology with a PDB object to check
        for consistency, just as gmx grompp does;
        if inconsistencies are found, prints a report
        :param pdb_object: a PDB instance # TODO either make our own or import from mdtraj, or make optional
        :return: int, 1 is consistent, 0 is not, -1 means naming incostistencies within the same element (C5T/C7 etc.)
        """
        pass
    
    def save_mod(self, outname):
        with open(outname, 'w') as outfile:
            for section in self.sections:
                for subsection in section.subsections:
                    outfile.write('\n[ {} ]\n'.format(subsection.write_header))
                    for line in subsection:
                        outfile.write(line)
