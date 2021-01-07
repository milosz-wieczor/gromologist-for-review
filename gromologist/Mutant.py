import gromologist as gml


class ProteinMutant:
    map_pro = {'ALA': 'A', 'CYS': 'C', 'CYX': 'C', 'CYM': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
               'HIS': 'H', 'HIE': 'H', 'HID': 'H', 'HSD': 'H', 'HSE': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
               'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
               'TRP': 'W', 'TYR': 'Y'}
    map_inv = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
               'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
               'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}

    def __init__(self, orig, target):
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
        self.mutate_from = []
        self.mutate_to = []
        length = max(len(orig_gps), len(target_gps))
        for i in range(length):
            orig = orig_gps[i] if i < len(orig_gps) else False
            targ = target_gps[i] if i < len(target_gps) else False
            if orig != targ:
                if orig and targ and self.compatible(orig, targ):
                    f, t = self.compatible_equivalents(orig, targ)
                    self.mutate_from.extend(f)
                    self.mutate_to.extend(t)
                if orig:
                    self.remove_from_orig.append(orig)
                if targ:
                    self.add_to_target.append(targ)

    def atoms_to_add(self):
        atoms = []
        hooks = []
        geo_refs = []
        bond_lengths = []
        topo_bonds = []
        for i in self.add_to_target:
            gatoms = self.groups(i)
            gbonds = self.bonds(i)
            ganchors = self.anchors(i)
            topo_bonds.extend(self.ring_closing_bonds(i))
            for n, j in enumerate(gatoms):
                if j not in self.mutate_to:
                    atoms.append(j)
                    geo_refs.append(ganchors[n])
                if gbonds[n][1] not in self.mutate_to:
                    hooks.append(gbonds[n][0])
                    bond_lengths.append(self.bond_lengths(gbonds[n]))
        return atoms, hooks, geo_refs, bond_lengths, topo_bonds

    def atoms_to_remove(self):
        atoms = []
        for i in self.remove_from_orig:
            gatoms = self.groups(i)
            for j in gatoms:
                if j not in self.mutate_from:
                    atoms.append(j)
        return atoms

    def atoms_to_mutate(self):
        return self.mutate_from, self.mutate_to

    @staticmethod
    def groups(key):
        gdict = {'CB': ['CB', 'HB1', 'HB2'], 'CG': ['CG', 'HG1', 'HG2'], 'CD': ['CD', 'HD1', 'HD2'],
                 'HA': ['HA2'], 'HB': ['HB3'], 'HD': ['HD3'], 'BM': ['CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23'],
                 'SH': ['SG', 'HG'], 'OH': ['OG', 'HG'], 'CO': ['CG', 'OD1', 'OD2'], 'DO': ['CD', 'OE1', 'OE2'],
                 'AM': ['CG', 'OD1', 'ND2', 'HD21', 'HD22'], 'AN': ['CD', 'OE1', 'NE2', 'HE21', 'HE22'],
                 'AR': ['CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ'], 'HG': ['HG3'],
                 'HZ': ['HZ'], 'LH': ['CG', 'HG', 'CD2', 'HD21', 'HD22', 'HD23'], 'OY': ['OH', 'HH'],
                 'SM': ['SD', 'CE', 'HE1', 'HE2', 'HE3'], 'KH': ['CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3'],
                 'HR': ['CG', 'ND1', 'CD2', 'HD2', 'CE1', 'NE2', 'HE1', 'HD1']
                 }
        return gdict[key]

    @staticmethod
    def bonds(key):
        bonds = {'CB': [('CA', 'CB'), ('CB', 'HB1'), ('CB', 'HB2')], 'CG': [('CB', 'CG'), ('CG', 'HG1'), ('CG', 'HG2')],
                 'CD': [('CG', 'CD'), ('CD', 'HD1'), ('CD', 'HD2')], 'BM': [('CA', 'CB'), ('CB', 'HB'), ('CB', 'CG2'),
                 ('CG2', 'HG21'), ('CG2', 'HG22'), ('CG2', 'HG23')], 'SH': [('CB', 'SG'), ('SG', 'HG')],
                 'OH': [('CB', 'OG'), ('OG', 'HG')], 'CO': [('CB', 'CG'), ('CG', 'OD1'), ('CG', 'OD2')],
                 'DO': [('CG', 'CD'), ('CD', 'OE1'), ('CD', 'OE2')], 'AM': [('CB', 'CG'), ('CG', 'OD1'), ('CG', 'ND2'),
                 ('ND2', 'HD21'), ('ND2', 'HD22')], 'AN': [('CG', 'CD'), ('CD', 'OE1'), ('CD', 'NE2'), ('NE2', 'HE21'),
                 ('NE2', 'HE22')], 'AR': [('CB', 'CG'), ('CG', 'CD1'), ('CD1', 'HD1'), ('CG', 'CD2'), ('CD2', 'HD2'),
                 ('CD1', 'CE1'), ('CE1', 'HE1'), ('CD2', 'CE2'), ('CE2', 'HE2'), ('CE2', 'CZ')], 'LH': [('CB', 'CG'),
                 ('CG', 'HG'), ('CG', 'CD2'), ('CD2', 'HD21'), ('CD2', 'HD22'), ('CD2', 'HD23')], 'OY': [('CZ', 'OH'),
                 ('OH', 'HH')], 'SM': [('CG', 'SD'), ('SD', 'CE'), ('CE', 'HE1'), ('CE', 'HE2'), ('CE', 'HE3')],
                 'KH': [('CD', 'CE'), ('CE', 'HE1'), ('CE', 'HE2'), ('CE', 'NZ'), ('NZ', 'HZ1'), ('NZ', 'HZ2'), ('NZ', 'HZ3')],
                 'HZ': [('CZ', 'HZ')], 'HA': [('CA', 'HA2')], 'HD': [('CD', 'HD3')], 'HG': [('CG', 'HG3')],
                 'HR': [('CB', 'CG'), ('CG', 'ND1'), ('CG', 'CD2'), ('CD2', 'HD2'), ('ND1', 'CE1'),
                        ('CD2', 'NE2'), ('CE1', 'HE1'), ('ND1', 'HD1')]}
        return bonds[key]

    @staticmethod
    def anchors(key):
        # obligatorily same order as in groups
        anchors = {'CB': [['CA', 'HA', 'C', 'N'], ['N', 'CA'], ['HA', 'CA']], 'CG': [['CB', 'CA', 'HB1', ('HB2', 'CG2')], [('HB2', 'CG2'), 'CB'],
                   ['HB1', 'CB']], 'CD': [['CG', 'CB', ('HG1', 'HG'), ('HG2', 'CD2')], [('HG2', 'CD2'), 'CG'], [('HG1', 'HG'), 'CG']],
                   'HA': [['CA', 'N', 'C', ('HA1', 'HA')]], 'HB': [['CB', 'CA', 'HB1', 'HB2']], 'HD': [['CD', 'CG', 'HD1', 'HD2']],
                   'BM': [['CA', 'HA', 'C', 'N'], ['N', 'CA'], ['HA', 'CA'], ['CA', 'CB'], ['CA', 'C'], ['CA', 'N']],
                   'SH': [['CB', 'CA', 'HB1', 'HB2'], ['CA', 'CB']], 'OH': [['CB', 'CA', 'HB1', ('HB2', 'CG2')], ['CA', 'CB']],
                   'CO': [['CB', 'CA', 'HB1', 'HB2'], ['CA', 'CB'], ['CG', 'CG', 'CB', 'OD1']],
                   'DO': [['CG', 'CB', 'HG1', 'HG2'], ['CB', 'CG'], ['CD', 'CD', 'CG', 'OE1']],
                   'AM': [['CB', 'CA', 'HB1', 'HB2'], ['CA', 'CB'], ['CG', 'CB', 'OD1', 'CG'], ['CB', 'CG'], ['OD1', 'CG']],
                   'AN': [['CG', 'CB', 'HG1', 'HG2'], ['CB', 'CG'], ['CD', 'CG', 'OE1', 'CD'], ['CG', 'CD'], ['OE1', 'CD']],
                   'AR': [['CB', 'CA', 'HB1', 'HB2'], ['CB', 'HB1', 'HB1', 'CA', 'CG', 'CG', 'HB1', 'HB1', 'CA', 'CG'], ['CG', 'CD1', 'CB'],
                          ['CB', 'HB2', 'HB2', 'CA', 'CG', 'CG', 'HB2', 'HB2', 'CA', 'CG'], ['CG', 'CD2', 'CB'], ['CB', 'CG'], ['CG', 'CD1'],
                          ['CB', 'CG'], ['CG', 'CD2'], ['CG', 'CD1']],
                   'HG': [['CG', 'CB', 'HG1', 'HG2']], 'OY': [['CG', 'CZ'], ['CA', 'CB']],
                   'HZ': [['CG', 'CZ']],
                   'LH': [['CB', 'CA', 'HB1', 'HB2'], ['HB2', 'CB'], ['HB1', 'CB'], ['CB', 'CG'], ['CB', 'CA'], ['CB', 'HB2']],
                   'SM': [['CA', 'CB'], ['CB', 'CG'], ['CG', 'SD'], ['CG', 'HG1'], ['CG', 'HG2']],
                   'KH': [['CB', 'CG'], ['CG', 'HG1'], ['CG', 'HG2'], ['CG', 'CD'], ['CD', 'CE'], ['CD', 'HD1'], ['CD', 'HD2']],
                   'HR': [['CB', 'CA', 'HB1', 'HB2'], ['CB', 'HB1', 'HB1', 'CA', 'CG', 'CG', 'HB1', 'HB1', 'CA', 'CG'],
                          ['CB', 'HB2', 'HB2', 'CA', 'CG', 'CG', 'HB2', 'HB2', 'CA', 'CG'], ['CG', 'CD2', 'CB'],
                          ['CB', 'CG', 'CD2'], ['CB', 'CG', 'ND1'], ['CE1', 'CE1', 'ND1', 'NE2'], ['ND1', 'ND1', 'CG', 'CE1']]}
        return anchors[key]

    @staticmethod
    def ring_closing_bonds(key):
        bonds = {'AR': [('CE2', 'CZ')], 'HR': [('CE1', 'NE2')]}
        return bonds[key]

    @staticmethod
    def aminoacids(key):
        aa = {'A': ['CB', 'HB'], 'C': ['CB', 'SH'], 'D': ['CB', 'CO'], 'E': ['CB', 'CG', 'DO'], 'F': ['CB', 'AR', 'HZ'],
              'G': ['HA'], 'H': ['CB', 'HR'], 'I': ['BM', 'CG', 'CD', 'HD'], 'K': ['CB', 'CG', 'CD', 'KH'],
              'L': ['CB', 'LH', 'CD', 'HD'], 'M': ['CB', 'CG', 'SM'], 'N': ['CB', 'AM'], 'Q': ['CB', 'CG', 'AN'],
              'R': ['CB', 'CG', 'CD', 'RH'], 'S': ['CB', 'OH'], 'T': ['BM', 'OH'], 'V': ['BM', 'CG', 'HG'],
              'Y': ['CB', 'AR', 'OY']}
        return aa[key]

    @staticmethod
    def compatible(k1, k2):
        comps = {('CB', 'BM'), ('CG', 'SH'), ('CG', 'OH'), ('CG', 'LH'), ('CO', 'AM'), ('AM', 'SH'), ('AM', 'OH'),
                 ('CO', 'SH'), ('CO', 'OH'), ('AN', 'DO'), }
        if (k1, k2) in comps or (k2, k1) in comps:
            return True
        return False

    @staticmethod
    def compatible_equivalents(k1, k2):
        ar = {('CB', 'BM'): [('CB', 'CB'), ('HB1', 'HB')], ('CG', 'SH'): [('CG', 'SG'), ('HG1', 'HG')],
              ('CG', 'OH'): [('CG', 'OG'), ('HG1', 'HG')], ('CG', 'LH'): [('CG', 'CG'), ('HG1', 'HG')],
              ('CO', 'AM'): [('CG', 'CG'), ('OD1', 'OD1'), ('OD2', 'ND2')],
              ('DO', 'AN'): [('CD', 'CD'), ('OE1', 'OE1'), ('OE2', 'NE2')], ('AM', 'SH'): [('CG', 'SG')],
              ('AM', 'OH'): [('CG', 'OG')], ('CO', 'SH'): [('CG', 'SG')], ('CO', 'OH'): [('CG', 'OG')], }
        try:
            x = ar[(k1, k2)]
            return [a[0] for a in x], [a[1] for a in x]
        except KeyError:
            x = ar[(k2, k1)]
            return [a[1] for a in x], [a[0] for a in x]

    @staticmethod
    def bond_lengths(bonded_pair):
        e1 = bonded_pair[0][0]
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
