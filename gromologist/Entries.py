

class Entry:
    """
    A generic class representing a single line in the topology.
    In an entry, the actual content and the comments are kept
    in two separate variables.
    """
    def __init__(self, content, subsection):
        self.subsection = subsection
        semicol_index = content.find(';')
        if semicol_index >= 0:
            self.content = content[:semicol_index].strip().split()
            self.comment = content[semicol_index:]
        else:
            self.content = content.strip().split()
            self.comment = ''
    
    @staticmethod
    def float_fmt(flt, fields=11, dpmax=8):
        """
        When a float of unknown precision is read, we do not want
        to clip off significant digits, but neither do we want
        to add too many decimal points. This function calculates
        how many dps we need to keep not to lose precision when
        handling ff params (by default, we clip to 8).
        :param flt: float, the number to be formatted
        :param fields: how many fields do we need overall in the fmt specifier
        :return: str, format specifier
        """
        nf = len(str(flt).split('.')[1])
        if nf > dpmax:
            nf = dpmax
        return "{:>" + str(fields) + "." + str(nf) + "f}"
        
    def __bool__(self):
        if not self.content and not self.comment:
            return False
        return True
    
    def __getitem__(self, item):
        return self.content[item]
    
    def __str__(self):
        """
        Fallback if no explicit formatting is implemented
        :return: str
        """
        return ' '.join(self.content) + self.comment + '\n'
    
    
class EntryBonded(Entry):
    """
    This Entry subclass is intended for entries that correspond
    to bonded interaction (bonds, pairs, angles, dihedrals)
    between specific atoms in the topology
    """
    fstr_suff = {('bonds', '1'): (float, float),
                 ('angles', '1'): (float, float),
                 ('angles', '5'): (float, float, float, float),
                 ('dihedrals', '9'): (float, float, int),
                 ('dihedrals', '4'): (float, float, int),
                 ('dihedrals', '2'): (float, float),
                 ('cmap', '1'): (float, float, float, float, float),
                 ('position_restraints', '1'): (int, int, int)}
    
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        self.atoms_per_entry = type(self.subsection).n_atoms[self.subsection.header]
        self.atom_numbers = tuple([int(x) for x in self.content[:self.atoms_per_entry]])
        self.interaction_type = self.content[self.atoms_per_entry]
        # type assignment should only be performed when asked to, i.e. outside of constructor, with read_types
        self.types_state_a = None
        self.types_state_b = None
        self.params_state_a = [] if len(self.content) == self.atoms_per_entry + 1 \
            else self.parse_params_state_a(self.content[self.atoms_per_entry + 1:])
        self.params_state_b = []  # TODO checks for int_type and self.subsection.header to determine if needed
        self.fstring = " ".join("{:>5d}" for n in range(self.atoms_per_entry)) + " {:>5s}"
    
    def read_types(self):
        atoms_sub = self.subsection.section.get_subsection('atoms')
        atoms_sub._get_dicts()
        num_to_type_a = atoms_sub.num_to_type
        num_to_type_b = atoms_sub.num_to_type_b
        self.types_state_a = tuple(num_to_type_a[num] for num in self.atom_numbers)
        self.types_state_b = tuple(num_to_type_b[num] for num in self.atom_numbers)
    
    def parse_params_state_a(self, excess_params):
        try:
            types = EntryBonded.fstr_suff[(self.subsection.header, self.interaction_type)]
        except KeyError:
            return [float(x) for x in excess_params]
        else:
            excess_params = excess_params[:len(types)]
            try:
                return [t(p) for p, t in zip(excess_params, types)]
            except ValueError:
                try:
                    params = self.subsection.section.top.defines[excess_params[0]]
                except KeyError:
                    raise RuntimeError("Cannot process: ", excess_params)
                else:
                    return [float(x) for x in params]
    
    def __str__(self):
        if self.params_state_a:
            fmt_suff = ""
            for parm in self.params_state_a:
                if isinstance(parm, int):
                    fmt_suff = fmt_suff + "{:>6d} "
                elif isinstance(parm, float):
                    fmt_suff = fmt_suff + self.float_fmt(parm)
            fstring = self.fstring + fmt_suff
            return fstring.format(*self.atom_numbers, self.interaction_type, *self.params_state_a) + self.comment
        else:
            return self.fstring.format(*self.atom_numbers, self.interaction_type) + self.comment
        
        
class EntryParam(Entry):
    """
    This Entry subclass represents a line containing force field
    parameters, e.g. bondtypes, angletypes, cmaptypes, pairtypes etc.
    that map a set of atom types to a set of FF-specific values
    """
    def __init__(self, content, subsection, processed=False):
        super().__init__(content, subsection)
        self.atoms_per_entry = type(self.subsection).n_atoms[self.subsection.header]
        self.types = tuple(self.content[:self.atoms_per_entry])
        if self.subsection.header == 'cmaptypes' and processed:
            self.modifiers = self.content[self.atoms_per_entry + 1:self.atoms_per_entry + 3]
            self.params = [float(x) for x in self.content[self.atoms_per_entry + 3:]]
            self.interaction_type = self.content[self.atoms_per_entry]
        elif self.subsection.header == 'cmaptypes' and not processed:
            self.modifiers = []
            self.params = self.content[self.atoms_per_entry + 1:]
            self.interaction_type = self.content[self.atoms_per_entry]
        elif self.subsection.header == 'defaults':
            self.modifiers = self.content
            self.params = []
            self.interaction_type = ''
        elif self.subsection.header == 'atomtypes':
            self.modifiers = self.content[self.atoms_per_entry:self.atoms_per_entry + 4]
            self.params = [float(x) for x in self.content[self.atoms_per_entry + 4:]]
            self.interaction_type = ''
        else:
            self.params = [float(x) for x in self.content[self.atoms_per_entry + 1:]]
            self.modifiers = []
            self.interaction_type = self.content[self.atoms_per_entry]
        if self.subsection.header == 'dihedraltypes' and self.interaction_type in ('9', '4'):
            self.params[-1] = int(self.params[-1])
            
    def format(self):
        fmt = {('bondtypes', '1'): "{:>8s} {:>8s}{:>6s}{:>13.8f}{:>13.2f}",
               ('angletypes', '5'): "{:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}{:>13.8f}{:>13.2f}",
               ('angletypes', '1'): "{:>8s} {:>8s} {:>8s}{:>6s}{:>13.8f}{:>13.2f}",
               ('dihedraltypes', '9'): "{:>8s} {:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}{:>6d}",
               ('dihedraltypes', '4'): "{:>8s} {:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}{:>6d}",
               ('dihedraltypes', '2'): "{:>8s} {:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}",
               ('atomtypes', ''): "{:>6s}{}{:>6s}{:>13s}{:>9s}{:>3s}{:>16.12f}{:>9.5f}",
               ('pairtypes', '1'): "{:>8s} {:>8s}{:>3s}{:>16.12f}{:>16.12f}",
               ('nonbond_params', '1'): "{:>8s} {:>8s}{:>3s}{:>16.12f}{:>16.12f}",
               ('implicit_genborn_params', ''): " {:8s}{:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}"}
        if (self.subsection.header, self.interaction_type) in fmt.keys():
            return fmt[(self.subsection.header, self.interaction_type)]
        else:
            return None
    
    def match(self, ext_typelist, int_type):
        if len(ext_typelist) != len(self.types):
            return False
        if self.interaction_type == int_type:
            if ext_typelist[0] == self.types[0] or ext_typelist[1] == self.types[1]:
                if all(ext_typelist[i] == self.types[i] for i in range(len(self.types)) if self.types[i] !='X'):
                    return True
            elif ext_typelist[0] == self.types[-1] or ext_typelist[1] == self.types[-2]:
                if all(ext_typelist[i] == self.types[len(self.types)-i-1] for i in range(len(self.types))
                       if self.types[i] !='X'):
                    return True
        return False
    
    def __repr__(self):
        if len(self.params) <= 4:
            return "Parameters entry with atomtypes {}, interaction type {} " \
                   "and parameters {}".format(self.types,
                                              self.interaction_type,
                                              ', '.join([str(x) for x in self.params]))
        else:
            return "Parameters entry with atomtypes {}, interaction type {} " \
                   "and parameters {}...".format(self.types,
                                                 self.interaction_type,
                                                 ', '.join([str(x) for x in self.params[:4]]))
        
    def __str__(self):
        """
        For cmaptypes, we rearrange lines to retrieve the matrix
        format lost during read-in; for other entry types, we
        delegate formatting to Subsection.fmt
        :return:
        """
        if self.subsection.header == 'cmaptypes':
            first = ((8 * "{} ")[:-1] + "\\\n").format(*self.types, self.interaction_type, *self.modifiers)
            npar = len(self.params)
            last = '\\\n'.join([((10 * "{} ")[:-1]).format(*self.params[10*n:10*(n+1)]) for n in range(int(npar/10))])
            if 10 * int(npar/10) != npar:
                last = last + '\\\n' + \
                       (((npar-10*int(npar/10)) * "{} ")[:-1]).format(*self.params[10*int(npar/10):]) + '\n'
            return first + last
        elif self.format():
            return self.format().format(*self.types, self.interaction_type, *self.modifiers, *self.params) +\
                   self.comment + '\n'
        else:
            return super().__str__()

        
class EntryAtom(Entry):
    """
    This Entry subclass corresponds to atoms defined in
    the [ atoms ] section of each molecule
    """
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        try:
            self.num, self.type, self.resid, self.resname, self.atomname, _, self.charge, self.mass = self.content[:8]
        except ValueError:
            self.num, self.type, self.resid, self.resname, self.atomname, _, self.charge = self.content[:7]
            self.mass = 0  # TODO get mass as the default from the ffnonbonded 'atomtypes' section
        self.type_b, self.charge_b, self.mass_b = None, None, None
        self.num, self.resid = int(self.num), int(self.resid)
        self.charge, self.mass = float(self.charge), float(self.mass)
        len(str(self.charge).split('.')[1])
        self.fstring = "{:>6d}{:>11s}{:>7d}{:>7s}{:>7s}{:>7d}"
    
    def __str__(self):
        fstring = self.fstring + self.float_fmt(self.charge) + self.float_fmt(self.mass) + '   '
        return fstring.format(self.num, self.type, self.resid, self.resname, self.atomname, self.num,
                              self.charge, self.mass) + self.comment
