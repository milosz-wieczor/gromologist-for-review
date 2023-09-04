import gromologist as gml


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
            self.comment = ' ' + content[semicol_index:].rstrip()
        else:
            self.content = content.strip().split()
            self.comment = ''
    
    @staticmethod
    def float_fmt(flt: float, fields: int = 11, dpmax: int = 8):
        """
        When a float of unknown precision is read, we do not want
        to clip off significant digits, but neither do we want
        to add too many decimal places. This function calculates
        how many dps we need to keep not to lose precision when
        handling ff params (by default, we clip to 8).
        :param flt: float, the number to be formatted
        :param fields: how many fields do we need overall in the fmt specifier
        :param dpmax: default limit on the number of decimal places
        :return: str, format specifier
        """
        try:
            nf = len(str(flt).split('.')[1])
        except IndexError:
            nf = 3
        if nf > 12:
            flt = float(str(flt).split('.')[0] + '.' + str(flt).split('.')[1][:12])
            nf = len(str(flt).split('.')[1])
        if nf > dpmax:
            dpmax = nf
            fields = dpmax + 3 + len(str(int(flt)))
        else:
            dpmax = nf
        return "{:>" + str(fields) + "." + str(dpmax) + "f}"

    @property
    def is_comment(self) -> bool:
        if not self.content and self.comment:
            return True
        else:
            return False

    @staticmethod
    def infer_type(val) -> type:
        try:
            _ = int(val)
        except ValueError:
            try:
                _ = float(val)
            except ValueError:
                return str
            else:
                return float
        else:
            return int
        
    def __bool__(self) -> bool:
        if not self.content and not self.comment:
            return False
        return True
    
    def __getitem__(self, item):
        return self.content[item]

    def is_header(self) -> bool:
        if len(self.content) == 0:
            return False
        return True if (self.content[0].strip().startswith('[') and self.content[-1].strip().endswith(']')) else False
    
    def __str__(self) -> str:
        """
        Fallback if no explicit formatting is implemented
        :return: str
        """
        return ' '.join(self.content) + ' ' + self.comment + '\n'
    
    
class EntryBonded(Entry):
    """
    This Entry subclass is intended for entries that correspond
    to bonded interaction (bonds, pairs, angles, dihedrals)
    between specific atoms in the topology
    """
    fstr_suff = {('bonds', '1'): (float, float),
                 ('bonds', '2'): (float, float),
                 ('angles', '1'): (float, float),
                 ('angles', '2'): (float, float),
                 ('angles', '5'): (float, float, float, float),
                 ('angles', '10'): (float, float),
                 ('dihedrals', '9'): (float, float, int),
                 ('dihedrals', '4'): (float, float, int),
                 ('dihedrals', '1'): (float, float, int),
                 ('dihedrals', '3'): (float, float, float, float, float, float),
                 ('dihedrals', '2'): (float, float),
                 ('cmap', '1'): (float,),
                 ('position_restraints', '1'): (float, float, float),
                 ('dihedral_restraints', '1'): (float, float, float),
                 ('constraints', '2'): (float,),
                 ('virtual_sites2', '1'): (float,),
                 ('settles', '1'): (float, float)}

    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        if subsection.header == 'exclusions':
            self.atoms_per_entry = len(self.content) - 1
        else:
            self.atoms_per_entry = type(self.subsection).n_atoms[self.subsection.header]
        self.atom_numbers = tuple([int(x) for x in self.content[:self.atoms_per_entry]])
        self.interaction_type = self.content[self.atoms_per_entry]
        try:
            self.params_per_entry = len(EntryBonded.fstr_suff[(subsection.header, str(self.interaction_type))])
        except KeyError:
            self.params_per_entry = 0
            # type assignment should only be performed when asked to, i.e. outside of constructor, with read_types
        self.types_state_a = None
        self.types_state_b = None
        self.atom_names = None
        self.params_state_a = []
        self.params_state_a_entry = []
        self.params_state_b = []
        self.params_state_b_entry = []
        self.fstr_mod = []
        if len(self.content) > self.atoms_per_entry + 1:
            try:
                self.parse_bonded_params(self.content[self.atoms_per_entry + 1:])
            except Exception as e:
                print("While trying to process line {}, subsection {}:".format(content, self.subsection.header))
                raise e
        self.fstring = " ".join("{:>5d}" for _ in range(self.atoms_per_entry)) + " {:>5s}"

    def _fstr_suff(self, query):
        if self.fstr_mod:
            return self.fstr_mod
        else:
            return EntryBonded.fstr_suff[query]

    def explicit_defines(self):
        if self.params_state_a and isinstance(self.params_state_a[0], str):
            try:
                self.params_state_a = self.subsection.section.top.defines[self.params_state_a[0]]
                self.fstr_mod = [self.infer_type(x) for x in self.params_state_a]
            except:
                pass

    def read_types(self):
        atoms_sub = self.subsection.section.get_subsection('atoms')
        atoms_sub.get_dicts()
        num_to_type_a = atoms_sub.num_to_type
        num_to_type_b = atoms_sub.num_to_type_b
        num_to_name = atoms_sub.num_to_name
        self.types_state_a = tuple(num_to_type_a[num] for num in self.atom_numbers)
        types_state_b = tuple(num_to_type_b[num] for num in self.atom_numbers)
        self.types_state_b = types_state_b if types_state_b != self.types_state_a else None
        self.atom_names = tuple(num_to_name[num] for num in self.atom_numbers)

    def parse_bonded_params(self, excess_params):
        try:
            _ = EntryBonded.fstr_suff[(self.subsection.header, self.interaction_type)]
        except KeyError:
            print((self.subsection.header, self.interaction_type))
            raise RuntimeError("Line '{}' contains unrecognized parameters".format(self.content))
        else:
            # if len(excess_params) == 1 and len(types) > 1:
            #     try:
            #         params = self.subsection.section.top.defines[excess_params[0]]
            #     except KeyError:
            #         raise RuntimeError("Cannot process: ", excess_params)
            #     else:
            #         self.params_state_a = [types[n](prm) for n, prm in enumerate(params)]
            try:
                _ = [float(x) for x in excess_params]
            except ValueError:
                # self.fstr_mod = list(EntryBonded.fstr_suff[(self.subsection.header, self.interaction_type)])
                for n in range(len(excess_params)):
                    try:
                        _ = float(excess_params[n])
                    except ValueError:
                        if not excess_params[n] in self.subsection.section.top.defines.keys():
                            if not self.subsection.section.top.ignore_missing_defines:
                                raise RuntimeError(f'undefined parameter {excess_params[n]} was found, try setting a value '
                                                   f'by specifying define={{"{excess_params[n]}": value}} when '
                                                   f'initializing Top, or set ignore_missing_defines=True')
                            else:
                                self.subsection.section.top.defines.update({excess_params[n]: 0})
                self.fstr_mod = [self.infer_type(x) for x in excess_params]
            types = self._fstr_suff((self.subsection.header, self.interaction_type))
            if len(excess_params) == len(types):
                self.params_state_a = [types[n](prm) for n, prm in enumerate(excess_params[:len(types)])]
            elif len(excess_params) == 2 * len(types):
                self.params_state_a = [types[n](prm) for n, prm in enumerate(excess_params[:len(types)])]
                self.params_state_b = [types[n](prm) for n, prm in enumerate(excess_params[len(types):])]
            else:
                raise RuntimeError("Cannot process: ", excess_params)

    def __str__(self):
        fmt_suff = ""
        for params in [self.params_state_a, self.params_state_b]:
            for parm in params:
                if isinstance(parm, int):
                    fmt_suff = fmt_suff + "{:>6d} "
                elif isinstance(parm, float):
                    fmt_suff = ' ' + fmt_suff + self.float_fmt(parm) + ' '
                elif isinstance(parm, str):
                    if len(parm) > 14:
                        fmt_suff = fmt_suff + "{{:>{}s}} ".format(len(parm)+2)
                    else:
                        fmt_suff = fmt_suff + "{:>15s} "
        fstring = self.fstring + fmt_suff
        return fstring.format(*self.atom_numbers, self.interaction_type, *self.params_state_a, *self.params_state_b) \
            + ' ' + self.comment + '\n'

        
class EntryParam(Entry):
    """
    This Entry subclass represents a line containing force field
    parameters, e.g. bondtypes, angletypes, cmaptypes, pairtypes etc.
    that map a set of atom types to a set of FF-specific values
    """
    def __init__(self, content: str, subsection: "gml.Subsection", processed: bool = False):
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
        if self.subsection.header == 'dihedraltypes':
            if any([self.infer_type(x) == float for x in self.types]):
                self.types = self.content[:2]
                self.interaction_type = self.content[2]
                self.params = [float(x) for x in self.content[3:]]
        if self.subsection.header == 'dihedraltypes' and self.interaction_type in ('9', '4', '1'):
            self.params[-1] = int(self.params[-1])
        self.identifier = self.subsection.header + '-' + '-'.join(self.types) + '-' + self.interaction_type
            
    def format(self) -> str:
        fmt = {('bondtypes', '1'): "{:>8s} {:>8s} {:>6s} {:>13.8f} {:>13.2f} ",
               ('angletypes', '5'): "{:>8s} {:>8s} {:>8s} {:>6s} {:>13.6f} {:>13.6f} {:>13.8f} {:>13.2f} ",
               ('angletypes', '1'): "{:>8s} {:>8s} {:>8s} {:>6s} {:>13.8f} {:>13.2f} ",
               ('dihedraltypes', '9'): "{:>8s} {:>8s} {:>8s} {:>8s} {:>6s} {:>13.6f} {:>13.6f} {:>6d} ",
               ('dihedraltypes', '4'): "{:>8s} {:>8s} {:>8s} {:>8s} {:>6s} {:>13.6f} {:>13.6f} {:>6d} ",
               ('dihedraltypes', '3'): "{:>8s} {:>8s} {:>8s} {:>8s} {:>6s} {:>13.6f} {:>13.6f} {:>13.6f} {:>13.6f} "
                                       "{:>13.6f} {:>13.6f} ",
               ('dihedraltypes', '2'): "{:>8s} {:>8s} {:>8s} {:>8s} {:>6s} {:>13.6f} {:>13.6f} ",
               ('dihedraltypes', '1'): "{:>8s} {:>8s} {:>8s} {:>8s} {:>6s} {:>13.6f} {:>13.6f} {:>6d} ",
               ('atomtypes', ''): "{:>6s} {} {:>6s} {:>13s} {:>9s} {:>3s} {:>16.12f} {:>9.5f} ",
               ('pairtypes', '1'): "{:>8s} {:>8s} {:>3s} {:>16.12f} {:>16.12f} ",
               ('nonbond_params', '1'): "{:>8s} {:>8s} {:>3s} {:>20.16f} {:>20.16f} ",
               ('implicit_genborn_params', ''): " {:8s} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} "}
        if (self.subsection.header, self.interaction_type) in fmt.keys():
            return fmt[(self.subsection.header, self.interaction_type)]
        else:
            return ''
    
    def match(self, ext_typelist: list, int_type: str) -> bool:
        """
        Checks if the entry matches an external atomtype list, given the interaction type
        :param ext_typelist: list of str, types to be checked
        :param int_type: str, integer identifier of the interaction type
        :return: bool, whether the entry matches the query
        """
        if not ext_typelist or len(ext_typelist) != len(self.types):
            return False
        if self.interaction_type == int_type:
            if ext_typelist[0] == self.types[0] or ext_typelist[1] == self.types[1] \
                    or ext_typelist[-1] == self.types[-1] or ext_typelist[-2] == self.types[-2]:
                if all(ext_typelist[i] == self.types[i] for i in range(len(self.types)) if self.types[i] != 'X'):
                    return True
            if ext_typelist[0] == self.types[-1] or ext_typelist[1] == self.types[-2] \
                    or ext_typelist[-1] == self.types[0] or ext_typelist[-2] == self.types[1]:
                if all(ext_typelist[i] == self.types[len(self.types)-i-1] for i in range(len(self.types))
                       if self.types[i] != 'X'):
                    return True
        return False
    
    def __repr__(self) -> str:
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
        
    def __str__(self) -> str:
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
            try:
                return self.format().format(*self.types, self.interaction_type, *self.modifiers, *self.params) +\
                       ' ' + self.comment + '\n'
            except:
                print((*self.types, self.interaction_type, *self.modifiers, *self.params))
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
            atomtypes = self.subsection.section.top.parameters.get_subsection('atomtypes')
            matching = [atype for atype in atomtypes if isinstance(atype, EntryParam) and atype.types[0] == self.type]
            try:
                self.mass = float(matching[0].modifiers[1])
            except IndexError:
                self.mass = 0
        self.num, self.resid = int(self.num), int(self.resid)
        self.charge, self.mass = float(self.charge), float(self.mass)
        if len(self.content) == 11:
            self.type_b, self.charge_b, self.mass_b = self.content[8], float(self.content[9]), float(self.content[10])
        else:
            self.type_b, self.charge_b, self.mass_b = None, None, None
        self.fstring = "{:>6d} {:>11s} {:>7d}{:>7s}{:>7s}{:>7d}"

    def __eq__(self, other):
        if not isinstance(other, gml.EntryAtom):
            return False
        for attr in ['atomname', 'num', 'resid', 'resname', 'charge', 'mass', 'type', 'type_b', 'charge_b', 'mass_b', 'molname']:
            if self.__getattribute__(attr) != other.__getattribute__(attr):
                return False
        return True

    @property
    def molname(self) -> str:
        """
        returns the name of the molecule this atom is part of
        :return: str, name of the molecule
        """
        return self.subsection.section.mol_name

    @property
    def element(self) -> str:
        return [x for x in self.atomname if not x.isdigit()][0]

    @property
    def ish(self, refstate: str = 'A') -> bool:
        """
        Tells if an atom is a hydrogen (useful e.g. for hydrogen mass repartitioning)
        WARNING not very relevant but might make mistakes with helium
        :param refstate: 'A' or 'B', alchemical state to take into account
        :return: bool, whether the atom is a hydrogen
        """
        if refstate == 'A':
            typecheck = self.type[0].upper() == 'H'
        elif refstate == 'B':
            typecheck = self.type_b[0].upper() == 'H'
        else:
            raise RuntimeError("refstate should be 'A' or 'B'")
        return typecheck

    @property
    def sigma(self) -> float:
        """
        Returns the LJ sigma value of an atom (in nm)
        :return: float, the value of sigma
        """
        entry = self._get_atomtype_entry()
        return float(entry.params[0])

    @property
    def epsilon(self) -> float:
        """
        Returns the LJ epsilon value of an atom (in kJ/mol)
        :return: float, the value of epsilon
        """
        entry = self._get_atomtype_entry()
        return float(entry.params[1])

    def _get_atomtype_entry(self):
        atomtypes = self.subsection.section.top.parameters.get_subsection('atomtypes')
        try:
            return [e for e in atomtypes if isinstance(e, EntryParam) and e.types[0] == self.type][-1]
        except IndexError:
            raise RuntimeError(f"Couldn't find non-bonded parameters for atomtype {self.type}")

    def __str__(self):
        fstring = self.fstring + self.float_fmt(self.charge) + self.float_fmt(self.mass) + '   '
        if self.type_b:
            alch_fstring = "{:>11s}" + self.float_fmt(self.charge_b) + self.float_fmt(self.mass_b)
            fstring += alch_fstring
            return fstring.format(self.num, self.type, self.resid, self.resname, self.atomname, self.num,
                                  self.charge, self.mass, self.type_b, self.charge_b, self.mass_b) + self.comment + '\n'
        else:
            return fstring.format(self.num, self.type, self.resid, self.resname, self.atomname, self.num,
                                  self.charge, self.mass) + ' ' + self.comment + '\n'
