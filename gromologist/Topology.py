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
        self.get_sections()
    
    def include_all(self):
        """
        includes all .itp files in the .top file to facilitate processing
        :return: None
        """
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
    
    def get_sections(self):
        special_sections = {'defaults', 'moleculetype', 'system'}
        special_lines = [n for n, l in enumerate(self.contents) if l.strip().strip('[]').strip() in special_sections]
        special_lines.append(len(self.contents) - 1)
        for beg, end in zip(special_lines[:-1], special_lines[1:]):
            self.sections.append(self.yield_sec(self.contents[beg:end]))
    
    def yield_sec(self, content):
        if 'defaults' in content[0]:
            return SectionParam(content, self)
        elif 'moleculetype' in content[0]:
            return SectionMol(content, self)
        else:
            return Section(content, self)
    
    def select_sections(self, section_name):
        return [s for s in self.sections if s.header == section_name]
    
    def parse_molecules(self):
        """
        Extracts individual molecules from topology, creating a dict that binds
        the moleculename to a list of numbers of sections that contain its bonded params
        :return: dict containing molname:[section_numbers] bindings
        """
        mols = {}
        # TODO this still needs some love
        # TODO need to rewrite so that mols collects SectionMol objects
        moltypes = self.select_sections('moleculetype')
        t, b, p, a, d = [self.select_sections(x) for x in ['atoms', 'bonds', 'pairs', 'angles', 'dihedrals']]
        bondeds = t + b + p + a + d
        moltypes.append(len(self.sections) - 1)
        for m, mnext in zip(moltypes[:-1], moltypes[1:]):
            sect = self.sections[m]
            molname = [x.split()[0] for x in sect if len(x.split()) > 1
                       and not x.startswith('[') and not x.startswith(';')][0]
            mols[molname] = [x for x in bondeds if m < x < mnext]
            mols[molname].sort()
        return mols
    
    def save_mod(self, outname):
        # TODO this also needs to be rewritten
        with open(outname, 'w') as outfile:
            for section in self.sections:
                for line in section:
                    outfile.write(line)
                outfile.write('\n')


