from typing import Tuple

class ProteinMutant:
    # W > Z -> not working
    map_pro = {'ALA': 'A', 'ASH': 'B', 'ASPP': 'B', 'CYS': 'C', 'CYX': 'C', 'CYM': 'C', 'ASP': 'D', 'GLU': 'E',
               'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'HIE': 'H', 'HID': 'H', 'HSD': 'H', 'HSE': 'H', 'ILE': 'I',
               'GLH': 'J', 'GLPP': 'J', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
               'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
    map_inv = {'A': 'ALA', 'B': 'ASH', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS',
               'I': 'ILE', 'J': 'GLH', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN',
               'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}

    def __init__(self, orig: str, target: str):
        try:
            target_1l = target if len(target) == 1 else ProteinMutant.map_pro[target]
        except KeyError:
            raise RuntimeError("Unrecognized residue specification: {}".format(target))
        try:
            self.name = ProteinMutant.map_pro[orig] + '2' + target_1l
        except KeyError:
            raise RuntimeError("You're trying to mutate residue {} which wasn't found in the database".format(orig))
        self.target_3l = ProteinMutant.map_inv[self.name[-1]]
        orig_gps = self.aminoacids(self.name[0])
        target_gps = self.aminoacids(self.name[-1])
        self.remove_from_orig = []
        self.add_to_target = []
        length = max(len(orig_gps), len(target_gps))
        diverged = False  # if one mismatch is found, side chain will be rebuilt all the way from there
        for i in range(length):
            orig = orig_gps[i] if i < len(orig_gps) else False
            targ = target_gps[i] if i < len(target_gps) else False
            if orig != targ or target_1l in 'IT' or diverged:
                if orig:
                    self.remove_from_orig.append(orig)
                if targ:
                    self.add_to_target.append(targ)
                diverged = True
        if self.name[-1] == 'P':
            self.remove_from_orig.append('PH')
        if self.name[0] == 'P':
            self.add_to_target.append('PH')

    def atoms_to_add(self) -> Tuple[list, list, list, list, list, list]:
        """
        For each atom group, identifies a set of atoms, bonds, anchors etc.
        that have to be included when implementing the mutation
        :return: tuple of lists, the required changes
        """
        # TODO make sure HA1/HA2 are properly assigned in GLY
        atoms = []
        hooks = []
        geo_refs = []
        bond_lengths = []
        topo_bonds = []
        afters = []
        for i in self.add_to_target:
            gatoms = self.groups(i)
            gbonds = self.bonds(i)
            ganchors = self.anchors(i)
            gafters = self.afters(i)
            topo_bonds.extend(self.ring_closing_bonds(i))
            for n, j in enumerate(gatoms):
                atoms.append(j)
                geo_refs.append(ganchors[n])
                hooks.append(gbonds[n][0])
                bond_lengths.append(self.bond_lengths(gbonds[n]))
                afters.append(gafters[n])
        return atoms, hooks, geo_refs, bond_lengths, topo_bonds, afters

    def atoms_to_remove(self) -> list:
        """
        For each atom group, identifies a set of atoms
        that have to be removed when implementing the mutation
        :return:
        """
        atoms = []
        for i in self.remove_from_orig:
            gatoms = self.groups(i)
            for j in gatoms:
                atoms.append(j)
        return atoms

    def check_chiral(self, struct, ca):
        """
        Makes sure the correct HA is removed in case of glycine
        :param struct: gml.Pdb instance
        :param ca: gml.Atom instance
        :return: None
        """
        if 'HA' in self.remove_from_orig:
            if not struct.check_chiral([ca], 'N', 'C', 'HA1', printing=False):
                resnum, resname, chain = ca.resnum, ca.resname, ca.chain
                chn = ' and chain ' + chain if chain.strip() else ''
                ha1 = struct.get_atom(f'name HA1 and resnum {resnum} and resname {resname} {chn}')
                ha2 = struct.get_atom(f'name HA2 and resnum {resnum} and resname {resname} {chn}')
                ha2.atomname = 'HA1'
                ha1.atomname = 'HA2'

    @staticmethod
    def groups(key: str) -> list:
        """
        Defines atoms constituting each atomic group
        (akin e.g. to a methyl or methylene group) in side chains
        :param key: str, name of the group
        :return: list of atoms in the group
        """
        gdict = {'CB': ['CB', 'HB1', 'HB2'],
                 'CG': ['CG', 'HG1', 'HG2'],
                 'CD': ['CD', 'HD1', 'HD2'],
                 'HA': ['HA2'],
                 'HB': ['HB3'],
                 'HD': ['HD3'],
                 'BM': ['CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23'],
                 'SH': ['SG', 'HG'],
                 'OH': ['OG', 'HG'],
                 'CO': ['CG', 'OD1', 'OD2'],
                 'DO': ['CD', 'OE1', 'OE2'],
                 'CP': ['CG', 'OD1', 'OD2', 'HD2'],
                 'DP': ['CD', 'OE1', 'OE2', 'HE2'],
                 'AM': ['CG', 'OD1', 'ND2', 'HD21', 'HD22'],
                 'AN': ['CD', 'OE1', 'NE2', 'HE21', 'HE22'],
                 'AR': ['CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ'],
                 'HG': ['HG3'],
                 'HZ': ['HZ'],
                 'LH': ['CG', 'HG', 'CD2', 'HD21', 'HD22', 'HD23'],
                 'OY': ['OH', 'HH'],
                 'SM': ['SD', 'CE', 'HE1', 'HE2', 'HE3'],
                 'KH': ['CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3'],
                 'HR': ['CG', 'ND1', 'CD2', 'HD2', 'CE1', 'NE2', 'HE1', 'HD1'],
                 'GG': ['CG1', 'HG11', 'HG12'],
                 'RH': ['NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22'],
                 'PD': ['CD', 'HD1', 'HD2'],
                 'PG': ['CG', 'HG1', 'HG2'],
                 'PB': ['CB', 'HB1', 'HB2'],
                 'PH': ['HN'],
                 'WR': ['CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CD2', 'CE3', 'HE3', 'CZ3', 'HZ3', 'CH2', 'HH2', 'CZ2',
                        'HZ2']
                 }
        return gdict[key]

    @staticmethod
    def afters(key: str) -> list:
        """
        For each atomic group, identifies by name atoms immediately
        preceding the ones to be introduced (allowing for variation)
        :param key: str, name of the group
        :return: list of atoms preceding the corresponding atom in self.groups[key] (or tuples if multiple names
        are accepted)
        """
        aftdict = {'CB': [('HA', 'HA2', 'HA1'), 'CB', 'HB1'],
                   'PB': [('HA', 'HA2', 'HA1'), 'CB', 'HB1'],
                   'CG': [('HB2', 'HB'), ('CG', 'CG1'), ('HG1', 'HG11')],
                   'PG': [('HB2', 'HB'), ('CG', 'CG1'), ('HG1', 'HG11')],
                   'CD': [('HG2', 'HG12', 'CD2'), ('CD', 'CD1'), ('HD1', 'HD11')],
                   'PD': [('HG2', 'HG12', 'CD2'), ('CD', 'CD1'), ('HD1', 'HD11')],
                   'HA': [('HA', 'HA1')],
                   'HB': ['HB2'],
                   'HD': [('HD2', 'HD12')],
                   'BM': ['HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22'],
                   'SH': ['HB2', 'SG'],
                   'OH': [('HB2', 'CG2'), ('OG', 'OG1')],
                   'CO': ['HB2', 'CG', 'OD1'],
                   'DO': [('HG2', ), 'CD', 'OE1'],
                   'CP': ['HB2', 'CG', 'OD1', 'OD2'],
                   'DP': [('HG2', ), 'CD', 'OE1', 'OE2'],
                   'AM': ['HB2', 'CG', 'OD1', 'ND2', 'HD21'],
                   'AN': ['HG2', 'CD', 'OE1', 'NE2', 'HE21'],
                   'GG': [('HB2', 'HG23'), 'CG1', 'HG11'],
                   'AR': ['HB2', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2'],
                   'HG': [('HG2', 'HG12')],
                   'HZ': ['CZ'], 'LH': ['HB2', 'CG', 'HG', 'CD2', 'HD21', 'HD22'],
                   'OY': ['CZ', 'OH'],
                   'SM': ['HG2', 'SD', 'CE', 'HE1', 'HE2'],
                   'KH': ['HD2', 'CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2'],
                   'HR': ['HB2', 'CG', 'ND1', 'CD2', 'HD2', 'CE1', 'NE2', 'HE1'],
                   'RH': ['HD2', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21'],
                   'PH': ['N'],
                   'WR': ['HB2', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CD2', 'CE3', 'HE3', 'CZ3', 'HZ3', 'CH2',
                          'HH2', 'CZ2']
                   }
        return aftdict[key]

    @staticmethod
    def bonds(key: str) -> list:
        """
        For each atomic group, identifies chemical bonds that connect atoms
        of that group, as well as the group to the previous one
        :param key: str, name of the group
        :return: list of lists, each containing a pair of atoms to be bonded (or tuples for ambiguous names)
        """
        bonds = {'CB': [['CA', 'CB'], ['CB', 'HB1'], ['CB', 'HB2']],
                 'PB': [['CA', 'CB'], ['CB', 'HB1'], ['CB', 'HB2']],
                 'CG': [['CB', 'CG'], [('CG', 'CG1'), 'HG1'], [('CG', 'CG1'), 'HG2']],
                 'PG': [['CB', 'CG'], [('CG', 'CG1'), 'HG1'], [('CG', 'CG1'), 'HG2']],
                 'GG': [['CB', 'CG1'], ['CG1', 'HG11'], ['CG1', 'HG12']],
                 'CD': [[('CG', 'CG1'), 'CD'], [('CD', 'CD1'), 'HD1'], [('CD', 'CD1'), 'HD2']],
                 'PD': [[('CG', 'CG1'), 'CD'], [('CD', 'CD1'), 'HD1'], [('CD', 'CD1'), 'HD2']],
                 'BM': [['CA', 'CB'], ['CB', 'HB'], ['CB', 'CG2'], ['CG2', 'HG21'], ['CG2', 'HG22'], ['CG2', 'HG23']],
                 'SH': [['CB', 'SG'], ['SG', 'HG']],
                 'OH': [['CB', 'OG'], [('OG', 'OG1'), 'HG']],
                 'CO': [['CB', 'CG'], ['CG', 'OD1'], ['CG', 'OD2']],
                 'DO': [['CG', 'CD'], ['CD', 'OE1'], ['CD', 'OE2']],
                 'CP': [['CB', 'CG'], ['CG', 'OD1'], ['CG', 'OD2'], ['OD2', 'HD2']],
                 'DP': [['CG', 'CD'], ['CD', 'OE1'], ['CD', 'OE2'], ['OE2', 'HE2']],
                 'AM': [['CB', 'CG'], ['CG', 'OD1'], ['CG', 'ND2'], ['ND2', 'HD21'], ['ND2', 'HD22']],
                 'AN': [['CG', 'CD'], ['CD', 'OE1'], ['CD', 'NE2'], ['NE2', 'HE21'], ['NE2', 'HE22']],
                 'AR': [['CB', 'CG'], ['CG', 'CD1'], ['CD1', 'HD1'], ['CG', 'CD2'], ['CD2', 'HD2'], ['CD1', 'CE1'],
                        ['CE1', 'HE1'], ['CD2', 'CE2'], ['CE2', 'HE2'], ['CE2', 'CZ']],
                 'LH': [['CB', 'CG'], ['CG', 'HG'], ['CG', 'CD2'], ['CD2', 'HD21'], ['CD2', 'HD22'], ['CD2', 'HD23']],
                 'OY': [['CZ', 'OH'], ['OH', 'HH']],
                 'SM': [['CG', 'SD'], ['SD', 'CE'], ['CE', 'HE1'], ['CE', 'HE2'], ['CE', 'HE3']],
                 'KH': [['CD', 'CE'], ['CE', 'HE1'], ['CE', 'HE2'], ['CE', 'NZ'], ['NZ', 'HZ1'], ['NZ', 'HZ2'],
                        ['NZ', 'HZ3']],
                 'PH': [['N', 'H']],
                 'HZ': [['CZ', 'HZ']],
                 'HA': [['CA', 'HA2']],
                 'HB': [['CB', 'HB3']],
                 'HD': [[('CD', 'CD1'), 'HD3']],
                 'HG': [[('CG', 'CG1'), 'HG3']],
                 'HR': [['CB', 'CG'], ['CG', 'ND1'], ['CG', 'CD2'], ['CD2', 'HD2'], ['ND1', 'CE1'], ['CD2', 'NE2'],
                        ['CE1', 'HE1'], ['ND1', 'HD1']],
                 'RH': [['CD', 'NE'], ['NE', 'HE'], ['NE', 'CZ'], ['CZ', 'NH1'], ['NH1', 'HH11'], ['NH1', 'HH12'],
                        ['CZ', 'NH2'], ['NH2', 'HH21'], ['NH2', 'HH22']],
                 'WR': [['CB', 'CG'], ['CG', 'CD1'], ['CD1', 'HD1'], ['CD1', 'NE1'], ['NE1', 'HE1'], ['NE1', 'CE2'],
                        ['CE2', 'CD2'], ['CD2', 'CE3'], ['CE3', 'HE3'], ['CE3', 'CZ3'], ['CZ3', 'HZ3'],
                        ['CZ3', 'CH2'], ['CH2', 'HH2'], ['CH2', 'CZ2'], ['CZ2', 'HZ2']]
                 }
        return bonds[key]

    @staticmethod
    def anchors(key: str) -> list:
        """
        For each atomic group, identifies sequences of atoms that define
        the geometry of the newly built group (through routines defined
        in gml.Pdb._vector)
        :param key: str, name of the group
        :return: list of lists
        """
        # obligatorily same order as in groups
        anchors = {'CB': [['CA', ('HA1', 'HA'), 'C', 'N'], ['N', 'CA'], [('HA1', 'HA'), 'CA']],
                   'PB': [['CA', ('HA1', 'HA'), 'C', 'N'], ['N', ('HA1', 'HA')],
                          [('HA1', 'HA'), 'CB', 'CB', 'CB', 'HB1', 'C', 'C', 'C']],
                   'CG': [['CB', 'CA', ('HB1', 'HB'), ('HB2', 'CG2')], [('HB2', 'CG2'), 'CB'], [('HB1', 'HB'), 'CB']],
                   'PG': [['CA', 'CB', 'CB', 'N', 'N', 'N'], [('HA1', 'HA'), 'CB'], ['C', 'CA']],
                   'GG': [['CB', 'CA', 'HB', 'CG2'], ['CG2', 'CB'], ['HB', 'CB']],
                   'CD': [[('CG', 'CG1'), 'CB', ('HG1', 'HG11', 'HG'), ('HG2', 'HG12', 'CD2')],
                          [('HG2', 'HG12', 'CD2'), ('CG', 'CG1')],
                          [('HG1', 'HG11', 'HG'), ('CG', 'CG1')]],
                   'PD': [['CB', 'N'], ['CA', 'N', 'C', 'CB', 'CA'], ['CD', 'N', 'CG', 'HD1']],
                   'HA': [['CA', 'N', 'C', ('HA1', 'HA')]],
                   'HB': [['CB', 'CA', 'HB1', 'HB2']],
                   'HD': [['CD', ('CG', 'CG1'), 'HD1', 'HD2']],
                   'BM': [['CA', ('HA1', 'HA'), 'C', 'N'], ['C', 'CA'], [('HA1', 'HA'), 'CA'], ['CA', 'CB'],
                          ['CA', 'C'], ['CA', 'N']],
                   'SH': [['CB', 'CA', 'HB1', 'HB2'], ['CA', 'CB']],
                   'OH': [['CB', 'CA', ('HB1', 'HB'), ('HB2', 'CG2')], ['CA', 'CB']],
                   'CO': [['CB', 'CA', 'HB1', 'HB2'], ['CA', 'CB'], ['CG', 'CG', 'CB', 'OD1']],
                   'DO': [['CG', 'CB', 'HG1', 'HG2'], ['CB', 'CG'], ['CD', 'CD', 'CG', 'OE1']],
                   'CP': [['CB', 'CA', 'HB1', 'HB2'], ['CA', 'CB'], ['CG', 'CG', 'CB', 'OD1'], ['CB', 'CG']],
                   'DP': [['CG', 'CB', 'HG1', 'HG2'], ['CB', 'CG'], ['CD', 'CD', 'CG', 'OE1'], ['CG', 'CD']],
                   'AM': [['CB', 'CA', 'HB1', 'HB2'], ['CA', 'CB'], ['CG', 'CB', 'OD1', 'CG'], ['CB', 'CG'],
                          ['OD1', 'CG']],
                   'AN': [['CG', 'CB', 'HG1', 'HG2'], ['CB', 'CG'], ['CD', 'CG', 'OE1', 'CD'], ['CG', 'CD'],
                          ['OE1', 'CD']],
                   'AR': [['CB', 'CA', 'HB1', 'HB2'], ['CB', 'HB1', 'HB1', 'CA', 'CG', 'CG', 'HB1', 'HB1', 'CA', 'CG'],
                          ['CG', 'CD1', 'CB'], ['CB', 'HB2', 'HB2', 'CA', 'CG', 'CG', 'HB2', 'HB2', 'CA', 'CG'],
                          ['CG', 'CD2', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CB', 'CG'], ['CG', 'CD2'], ['CG', 'CD1']],
                   'HG': [['CG', 'CB', ('HG1', 'HG11'), ('HG2', 'HG12')]],
                   'OY': [['CG', 'CZ'], ['CA', 'CB']],
                   'HZ': [['CG', 'CZ']],
                   'PH': [['CA', 'N', 'N', ('CB', 'HA2'), ('CB', 'HA2')]],
                   'LH': [['CB', 'CA', 'HB1', 'HB2'], ['HB2', 'CB'], ['HB1', 'CB'], ['CB', 'CG'], ['CB', 'CA'],
                          ['CB', 'HB2']],
                   'SM': [['CA', 'CB'], ['CB', 'CG'], ['CG', 'SD'], ['CG', 'HG1'], ['CG', 'HG2']],
                   'KH': [['CB', 'CG'], ['CG', 'HG1'], ['CG', 'HG2'], ['CG', 'CD'], ['CD', 'CE'], ['CD', 'HD1'],
                          ['CD', 'HD2']],
                   'HR': [['CB', 'CA', 'HB1', 'HB2'], ['CB', 'HB1', 'HB1', 'CA', 'CG', 'CG', 'HB1', 'HB1', 'CA', 'CG'],
                          ['CB', 'HB2', 'HB2', 'CA', 'CG', 'CG', 'HB2', 'HB2', 'CA', 'CG'], ['CG', 'CD2', 'CB'],
                          ['CB', 'CG', 'CD2'], ['CB', 'CG', 'ND1'], ['CE1', 'CE1', 'ND1', 'NE2'],
                          ['ND1', 'ND1', 'CG', 'CE1']],
                   'RH': [['CD', 'HD1', 'HD2', 'CG'], ['CD', 'NE', 'CG'], ['CG', 'CD'], ['CD', 'NE'], ['NE', 'CZ'],
                          ['NE', 'HE'], ['HE', 'NE'], ['NE', 'CZ'], ['NE', 'CD']],
                   # ['CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CD2', 'CE3', 'HE3', 'CZ3', 'HZ3', 'CZ2', 'HZ2',
                   # 'CH2', 'HH2']
                   'WR': [['CB', 'CA', 'HB1', 'HB2'],  # CG
                          ['CB', 'CG', 'CA'],  # CD1
                          ['CG', 'CB', 'CD1'],  # HD1
                          ['CD1', 'CG', 'HD1', 'HD1'],  # NE1
                          ['CG', 'CD1'],  # HE1
                          ['CG', 'CB', 'CD1', 'CD1'],  # CE2
                          ['HE1', 'CD1'],  # CD2
                          ['CD1', 'CD2'],  # CE3
                          ['CE2', 'CD2'],  # HE3
                          ['HD1', 'NE1'],  # CZ3
                          ['CE2', 'CZ3'],  # HZ3
                          ['CD2', 'CE2'],  # CH2
                          ['CD2', 'CH2'],  # HH2
                          ['CD2', 'CD1'],  # CZ2
                          ['CE3', 'CZ2']]  # HZ2
                   }
        return anchors[key]

    @staticmethod
    def ring_closing_bonds(key: str) -> list:
        """
        If a residue contains a ring, this fn defines an extra bond
        that closes the ring (might require special treatment)
        :param key: str, name of the group
        :return: list of tuples, each defines a ring-closing bond
        """
        bonds = {'AR': [('CE1', 'CZ')], 'HR': [('CE1', 'NE2')], 'PD': [('CD', 'N')],
                 'WR': [('CZ2', 'CE2'), ('CG', 'CD2')]}
        return bonds[key] if key in bonds.keys() else []

    @staticmethod
    def aminoacids(key: str) -> list:
        """
        Lists groups that make up individual amino acids
        :param key: str, 1-letter code of the amino acid
        :return: list, groups that define the side chain
        """
        aa = {'A': ['CB', 'HB'],
              'B': ['CB', 'CP'],  # protonated ASP
              'C': ['CB', 'SH'],
              'E': ['CB', 'CG', 'DO'],
              'D': ['CB', 'CO'],
              'F': ['CB', 'AR', 'HZ'],
              'G': ['HA'],
              'H': ['CB', 'HR'],
              'I': ['BM', 'GG', 'CD', 'HD'],
              'J': ['CB', 'CG', 'DP'],  # protonated GLU
              'K': ['CB', 'CG', 'CD', 'KH'],
              'L': ['CB', 'LH', 'CD', 'HD'],
              'M': ['CB', 'CG', 'SM'],
              'N': ['CB', 'AM'],
              'Q': ['CB', 'CG', 'AN'],
              'P': ['PB', 'PG', 'PD'],
              'R': ['CB', 'CG', 'CD', 'RH'],
              'S': ['CB', 'OH'],
              'T': ['BM', 'OH'],
              'V': ['BM', 'CG', 'HG'],
              'W': ['CB', 'WR'],
              'Y': ['CB', 'AR', 'OY']}
        return aa[key]

    @staticmethod
    def bond_lengths(bonded_pair: list) -> float:
        """
        Defines bond lengths of pairs of atoms
        :param bonded_pair: list of str, names of the elements that create the bond
        :return: float, length of the bond
        """
        e1 = bonded_pair[0][0] if not isinstance(bonded_pair[0], tuple) else bonded_pair[0][0][0]
        e2 = bonded_pair[1][0]
        if e1 in 'C' and e2 in 'C':
            return 1.5
        elif e1 in 'CO' and e2 in 'CO':
            return 1.4
        elif e1 in 'CN' and e2 in 'CN':
            return 1.45
        elif e1 in 'CH' and e2 in 'CH':
            return 1.1
        elif e1 in 'NOH' and e2 in 'NOH':
            return 1.0
        elif e1 in 'CS' and e2 in 'CS':
            return 1.8
        elif e1 in 'HS' and e2 in 'HS':
            return 1.35
        else:
            raise RuntimeError('Cannot assign bond length to pair: {} {}'.format(e1, e2))
