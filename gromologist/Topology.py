import os
from .Section import *
from collections import OrderedDict
from .Pdb import Pdb


class Top:
    def __init__(self, filename, gmx_dir='/usr/share/gromacs/top/', pdb=None, ignore_ifdef=False):
        """
        A class to represent and contain the Gromacs topology file and provide
        tools for editing topology elements
        :param filename: str, path to the .top file
        :param gmx_dir: str, Gromacs FF directory
        :param pdb: str, path to a matching PDB file
        :param ignore_ifdef: bool, whether to ignore #include statements within #ifdef blocks (e.g. posre.itp)
        """
        # TODO maybe allow for construction of a blank top with a possibility to read data later?
        self.pdb = None if pdb is None else Pdb(pdb, top=self)
        self.fname = filename
        self.top = self.fname.split('/')[-1]
        self.dir = os.getcwd() + '/' + '/'.join(self.fname.split('/')[:-1])
        self._contents = open(self.fname).readlines()
        self._gromacs_dir = gmx_dir
        self._include_all(ignore_ifdef)
        self.sections = []
        self._parse_sections()
        if self.top.endswith('top'):
            self.system, self.charge, self.natoms = self._read_system_properties()
        
    def __repr__(self):
        return "Topology with {} atoms and total charge {}".format(self.natoms, self.charge)
    
    def add_pdb(self, pdbfile):
        """
        Allows to pair a PDB file with the topology after the instance was initialized
        :param pdbfile: str, path to PDB file
        :return: None
        """
        self.pdb = Pdb(pdbfile, top=self)
    
    def _include_all(self, ign_ifdef):
        """
        includes all .itp files in the .top file to facilitate processing
        :return: None
        """
        # TODO optionally insert .mdp defines in constructor & pass here?
        ignore_lines = self._find_ifdef_lines() if ign_ifdef else set()
        lines = [i for i in range(len(self._contents)) if self._contents[i].strip().startswith("#include")
                 and i not in ignore_lines]
        while len(lines) > 0:
            lnum = lines[0]
            to_include, extra_prefix = self._find_in_path(self._contents.pop(lnum).split()[1].strip('"\''))
            contents = self._add_prefix_to_include(open(to_include).readlines(), extra_prefix)
            self._contents[lnum:lnum] = contents
            ignore_lines = self._find_ifdef_lines() if ign_ifdef else set()
            lines = [i for i in range(len(self._contents)) if self._contents[i].startswith("#include")
                     and i not in ignore_lines]
    
    def _find_ifdef_lines(self):
        """
        Finds #ifdef/#endif blocks in the topology if user
        explicitly asks to ignore them
        :return: set, int numbers of all lines that fall within an #ifdef/#endif block
        """
        ignore_set = set()
        counter = 0
        for n, line in enumerate(self._contents):
            if line.strip().startswith("#ifdef"):
                counter += 1
            elif line.strip().startswith("#endif"):
                counter -= 1
            if counter > 0:
                ignore_set.add(n)
        return ignore_set
    
    def clear_sections(self):
        """
        Removes all SectionMol instances that are not part
        of the system definition in [ system ]
        # TODO optionally we could also delete all unused params
        # TODO we could also have a fn that packs .top into .itps again based on section definitions
        :return: None
        """
        if self.system is None:
            raise AttributeError("System properties have not been read, this is likely not a complete .top file")
        sections_to_delete = []
        for section_num, section in enumerate(self.sections):
            if isinstance(section, SectionMol) and section.mol_name not in self.system.keys():
                sections_to_delete.append(section_num)
        self.sections = [s for n, s in enumerate(self.sections) if n not in sections_to_delete]
    
    def _find_in_path(self, filename):
        """
        looks for a file to be included in either the current directory
        or in Gromacs directories (as given by user), in order to
        include all .itp files in a single .top file
        :param filename: str, name of the file to be searched for
        :return: None
        """
        pref = ''
        if filename in os.listdir(self.dir):
            return self.dir + '/' + filename, pref
        else:
            if '/' in filename:
                pref = '/'.join(filename.split('/')[:-1])
            suff = filename.split('/')[-1]
            if pref.startswith('/') or pref.startswith('./'):
                first = pref.split('/')[1]
            else:
                first = pref.split('/')[0]
            if first in os.listdir(self.dir) and suff in os.listdir(self.dir + '/' + pref):
                return self.dir + '/' + pref + '/' + suff, pref
            elif suff in os.listdir(self._gromacs_dir + '/' + pref):
                return self._gromacs_dir + '/' + pref + '/' + suff, pref
        raise FileNotFoundError('file {} not found in neither local nor Gromacs directory.\n'
                                'If the file is included in an #ifdef block, please try setting'
                                ' ignore_ifdef=True'.format(filename))
    
    @staticmethod
    def _add_prefix_to_include(content, prefix):
        """
        Modifies #include statements if nested #includes
        point to different directories
        :param content: list of str, content of the included file
        :param prefix: str, directory name to add in the nested include
        :return: list of str, modified content
        """
        if prefix:
            for nline, line in enumerate(content):
                if line.strip().startswith("#include"):
                    try:
                        index = line.index('"')
                    except ValueError:
                        index = line.index("'")
                    newline = line[:index+1] + prefix + '/' + line[index+1:]
                    content[nline] = newline
        return content
    
    def _parse_sections(self):
        """
        Cuts the content in sections as defined by the position
        of special headers, and builds the self.sections list
        :return: None
        """
        special_sections = {'defaults', 'moleculetype', 'system'}
        special_lines = [n for n, l in enumerate(self._contents)
                         if l.strip() and l.strip().strip('[]').strip().split()[0] in special_sections]
        special_lines.append(len(self._contents))
        for beg, end in zip(special_lines[:-1], special_lines[1:]):
            self.sections.append(self._yield_sec(self._contents[beg:end]))
    
    def _yield_sec(self, content):
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
        
    def _read_system_properties(self):
        """
        Reads in system composition based on the [ molecules ] section
        and calculates the number of atoms and charge of the system
        :return: OrderedDict, contains molecule_name:number_of_molecules pairs and preserves the order read from file
                 float, total charge of the system
                 int, number of atoms in the system
        """
        system_subsection = [s.get_subsection('molecules') for s in self.sections
                             if 'molecules' in [ss.header for ss in s.subsections]]
        molecules = OrderedDict()  # we want to preserve the order of molecules in the system for e.g. PDB checking
        natoms, charge = 0, 0
        if len(system_subsection) == 0:
            raise KeyError
        elif len(system_subsection) > 1:
            print("Section 'molecules' not present in the topology")
            return None, None, None
        for e in system_subsection[0]:
            if e.content:
                molecules[e.content[0]] = int(e.content[1])
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
    
    def check_pdb(self):
        """
        Compares the topology with a PDB object to check
        for consistency, just as gmx grompp does;
        if inconsistencies are found, prints a report
        :return: None
        """
        if self.pdb:
            self.pdb.check_top()
        else:
            raise AttributeError("No PDB file has been bound to this topology")
    
    def save_top(self, outname):
        """
        Saves the combined topology to the specified file
        :param outname: str, file name for output
        :return: None
        """
        with open(outname, 'w') as outfile:
            for section in self.sections:
                for subsection in section.subsections:
                    outfile.write('\n[ {} ]\n'.format(subsection.write_header))
                    for line in subsection:
                        outfile.write(line)
