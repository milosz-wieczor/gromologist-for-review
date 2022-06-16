import os
import gromologist as gml


class Mod:
    counter = 0  # only used for naming with unique numbers

    def __init__(self, top, master, resnr='', resnames='', names=''):
        """
        Base class to work with parameter modifications.
        The class is actually used for the calculation of original
        (unmodified) quantities, and is subclassed by ModParam,
        ModNbfix and ModAtom.
        :param top: list of Top objects corresponding to the modified topology
        :param resnr: a list of tuples, each containing a list of 3 or 4 residue numbers
        :param resnames: a list of tuples, each containing a list of 3 or 4 residue names
        :param names: a list of tuples, each containing a list of 3 or 4 atom names
        """
        Mod.counter += 1
        self.td = master
        self.counter = Mod.counter
        if isinstance(top, str):
            self.top = [gml.Top(top)]
        elif isinstance(top, gml.Top):
            self.top = [top]
        elif isinstance(top[0], str):
            self.top = [gml.Top(t) for t in top]
        else:
            self.top = top
        assert isinstance(self.top[0], gml.Top)
        self.data = None
        self.diff = None
        self.topname = ''
        self.trajs = []
        self.dpar = 0.001
        self.sigma, self.chg, self.eps, self.nbs, self.nbe, self.dih, self.ang, self.nothing = [''] * 8
        assert bool(names) == bool(resnames) == bool(resnr)
        self.names = names
        self.res = resnr
        self.resnames = resnames
        assert len(self.names) == len(self.res) == len(self.resnames) > 0

    def __str__(self):
        """
        string representation of a Mod object, should be unique
        :return: mod one-letter code + Mod index, 1-based
        """
        return "{}-{}".format(''.join([self.sigma, self.chg, self.eps, self.nbs, self.nbe, self.dih, self.ang,
                                       self.nothing]), self.counter)

    def __eq__(self, other):
        """
        assumes unique numbering of instances
        :param other: another Mod (or subclass) instance to compare with
        :return: bool
        """
        return str(self) == str(other)

    def save_mod(self, path, name):
        """
        wraps Top's save_mod to catch the file name
        :param path: str, filename to save .top to
        :param name: str, name of the file
        :return: None
        """
        for t in self.top:
            self.topname = path
            final_name = self.topname + name + '-' + t.top
            if not os.path.exists(final_name):
                t.save_mod(final_name)
            else:
                print(f'file {final_name} already exists, leaving as it is')

    def goto_mydir(self):
        """
        shortcut to the ModAtom's working directory
        :return: None
        """
        os.chdir(self.topname)


class ModAtom(Mod):
    def __init__(self, top, master, resnr='', resnames='', names='', changes=''):
        """
        Subclass to work with modified atomic params (that is,
        sigma and epsilon defined for individual types, or charge defined
        for individual atoms).
        The class takes in a set of atoms that all correspond
        to a single atom type (has to be the same for all atoms);
        these atoms' type will be cloned from "type" to "Ytype",
        and the corresponding sigma/epsilon/charge parameters will be modified.
        :param top: a Top object corresponding to the modified topology
        :param resnr: a list of residue numbers
        :param resnames: a list of residue names
        :param names: a list of of atom names
        :param changes: single character denoting the change to be performed
        """
        super(ModAtom, self).__init__(top, master, resnr, resnames, names)
        self.sigma = 's' if 's' in changes else ''
        self.eps = 'e' if 'e' in changes else ''
        self.chg = 'c' if 'c' in changes else ''
        if not (self.chg or self.eps or self.sigma):
            raise ValueError("no mode selected appropriate for an ModAtom object")
        if names and resnr:
            for t in self.top:
                self.types_are_consistent(t)
                self.type = [m.get_atoms(f'atomname {self.names[0]} and resid {self.res[0]} and '
                                         f'resname {self.resnames[0]}').type for m in t.molecules][0][0]
                if not self.chg:
                    t.params.clone_type(atomtype=self.type, prefix='Y')
                self.mod_sigma_eps(t)
                self.mod_chg(t)
        else:
            self.names = ['orig']
            self.res = 0

    def mod_chg(self, top):
        """
        looks for a line in topology section 'atoms'
        and modifies it as needed
        :return: None
        """
        for x in range(len(self.names)):
            name = self.names[x]
            resname = self.resnames[x]
            resnum = self.res[x]
            for section in top.list_sections('atoms'):
                for l in range(len(top.sections[section])):
                    lspl = top.sections[section][l].split()
                    if len(lspl) >= 7 and lspl[4] == name and int(lspl[2]) == resnum and lspl[3] == resname:
                        modline = self.increment_charge_line(lspl)
                        top.sections[section][l] = modline

    def increment_charge_line(self, cont_list):
        """
        takes a line from the topology section 'atoms'
        and increases charge by self.dpar;
        also changes atomtype to modified (Y-prefixed)
        :param cont_list: line passed as a list of non-whitespace strings
        :return: an assembled line with type and/or charge modified
        """
        cont_list[6] = float(cont_list[6])
        if len(cont_list) == 7:
            cont_list = cont_list + [' ']
        else:
            cont_list = cont_list[:8]
        if self.chg:
            cont_list.extend([cont_list[1], float(cont_list[6]) + self.dpar, cont_list[7]])
        if self.sigma or self.eps:
            cont_list = cont_list[:8]
            cont_list.append('Y' + cont_list[1])
            cont_list.append(cont_list[6])
            cont_list.append(cont_list[7])
        fstring = '{:>6s}{:>11s}{:>7s}{:>7s}{:>7s}{:>7s}{:>11.4f}{:>11s}{:>7s}{:>11.4f}{:>11s}\n'
        return fstring.format(*cont_list[:11])

    def mod_sigma_eps(self, top):
        """
        looks for a line in topology section 'atomtypes'
        and modifies it as needed
        :return: None
        """
        section = top.list_sections('atomtypes')[0]
        for l in range(len(top.sections[section])):
            lspl = top.sections[section][l].split()
            if len(lspl) >= 6 and lspl[0] == 'Y' + self.type:
                modline = self.increment_sigma_eps_line(lspl)
                top.sections[section][l] = modline

    def increment_sigma_eps_line(self, cont_list):
        """
        takes a line from the topology section 'atomtypes'
        and increases sigma by self.dpar
        :param cont_list: line passed as a list of non-whitespace strings
        :return: an assembled line with sigma modified
        """
        pos_A = cont_list.index("A")
        if self.sigma:
            if pos_A == 4:
                cont_list[5] = float(cont_list[5]) + self.dpar
            elif pos_A == 3:
                cont_list[4] = float(cont_list[4]) + self.dpar
            else:
                raise RuntimeError("Improper line formating: {}".format(" ".join(cont_list)))
        if self.eps:
            if pos_A == 4:
                cont_list[6] = float(cont_list[6]) + self.dpar
            elif pos_A == 3:
                cont_list[5] = float(cont_list[5]) + self.dpar
            else:
                raise RuntimeError("Improper line formating: {}".format(" ".join(cont_list)))
        if pos_A == 4:
            fstring = '{:<11s}{:>3s}{:>11s}{:>10s}{:>3s}{:>20}{:>20}'
        elif pos_A == 3:
            fstring = '{:<11s}{:>11s}{:>10s}{:>3s}{:>20}{:>20}'
        return fstring.format(*cont_list[:7])

    def types_are_consistent(self, top):
        """
        if many atoms are passed in a single input line, need to make sure
        they are all of the same type
        :return: True if all are of the same type, False otherwise
        """
        consistent = all([self.get_type(top, self.names[0], self.res[0], self.resnames[0])
                          == self.get_type(top, self.names[n], self.res[n], self.resnames[n])
                          for n in range(len(self.names))])
        if not consistent:
            raise ValueError(
                "atoms within a single charge, sigma or epsilon modification need to have consistent types")


class ModNbfix(Mod):
    def __init__(self, top, master, resnr='', resnames='', names='', changes=''):
        """
        Subclass to work with modified NBFIX params (pairwise corrections,
        sigma and epsilon defined for pairs of types).
        The class allows two lists (actually, tuples) of atoms that correspond
        to two atom types (atomtypes have to be consistent within each tuple);
        these atoms' types will be cloned from (type1, type2) to (Xtype1,
        Ytype2), and the corresponding NBFIX parameters will be modified.
        The modification can be applied to all atoms bearing the specific type
        as well as only a subgroup of these.
        :param top: a Top object corresponding to the modified topology
        :param resnr: a list of two tuples, each containing a list of residue numbers
        :param resnames: a list of two tuples, each containing a list of residue names
        :param names: a list of two tuples, each containing a list of atom names
        :param changes: single character denoting the change to be performed
        """
        super(ModNbfix, self).__init__(top, master, resnr, resnames, names)
        self.nbs = 'n' if 'n' in changes else ''
        self.nbe = 'm' if 'm' in changes else ''
        if not (self.nbs or self.nbe):
            raise ValueError("no mode selected appropriate for an ModParam object")
        for t, tnames, tres, tresnames in zip(self.top, self.names, self.res, self.resnames):
            self.LJ_sigmas = []
            self.LJ_epsilons = []
            self.NB_sigmas = None
            self.NB_epsilons = None
            self.types_are_consistent(t, tnames, tres, tresnames)
            if 'nonbond_params' not in t.section_headers:
                self.add_nbfix_section(t)
            try:
                self.types = (self.get_type(t, tnames[0][0], tres[0][0], tresnames[0][0]),
                              self.get_type(t, tnames[1][0], tres[1][0], tresnames[1][0]))
            except IndexError:
                pass
            else:
                self.check_lj_params(t)
                if self.types[0] == self.types[1]:
                    self.prefixes = ['Y', 'Y']
                else:
                    self.prefixes = ['Q', 'Y']
                # we mod the first type
                self.clone_type(t, atomtype=self.types[0], prefix=self.prefixes[0])
                self.prefix_type(t, tnames, tres, tresnames, 0)
                # and then second, if different than first
                if self.types[0] != self.types[1]:
                    self.clone_type(t, atomtype=self.types[1], prefix=self.prefixes[1])
                    self.prefix_type(t, tnames, tres, tresnames, 1)
                self.sort_dihedrals(t)
                self.mod_nb_sigma_eps(t)

    def mod_nb_sigma_eps(self, top):
        """
        looks for a line in topology section 'nonbond_params'
        and modifies it as needed
        :return: None
        """
        for section in top.list_sections('nonbond_params'):
            to_remove = []
            for line in range(len(top.sections[section])):
                lspl = top.sections[section][line].split()
                prefixed_types = (self.prefixes[0] + self.types[0], self.prefixes[1] + self.types[1])
                if len(lspl) > 4 and ((lspl[0], lspl[1]) == prefixed_types or
                                      (lspl[1], lspl[0]) == prefixed_types):
                    to_remove.append(top.sections[section][line])
            for redundant in to_remove:
                top.sections[section].remove(redundant)
        section = top.list_sections('nonbond_params')[0]
        top.sections[section].extend(self.increment_nb_sigma_eps_line(self.NB_epsilons is None and
                                                                      self.NB_sigmas is None))
        #
        # for l in range(len(top.sections[section])):
        #     lspl = top.sections[section][l].split()
        #     if len(lspl) > 4 and ((lspl[0] == self.prefixes[0] + self.types[0]
        #                            and lspl[1] == self.prefixes[1] + self.types[1])
        #                           or (lspl[1] == self.prefixes[0] + self.types[0]
        #                               and lspl[0] == self.prefixes[1] + self.types[1])):
        #         modline = self.increment_nb_sigma_eps_line(lspl)
        #         top.sections[section][l] = modline

    def increment_nb_sigma_eps_line(self, add_original):
        """
        takes a line from the topology section 'atomtypes'
        and increases the specific value by self.dpar
        :return: an assembled line with sigma modified
        """
        fstring = '  {:<11s}{:<11s}{:>3s}{:>20.8f}{:>20.8f}\n'
        if self.NB_epsilons is None and self.NB_sigmas is None:
            sigma = 0.5 * (self.LJ_sigmas[0] + self.LJ_sigmas[1])
            epsilon = (self.LJ_epsilons[0] * self.LJ_epsilons[1]) ** 0.5
        else:
            sigma = self.NB_sigmas
            epsilon = self.NB_epsilons
        cont_list = [self.types[0], self.types[1], '1', sigma, epsilon]
        mod_list = [self.prefixes[0] + self.types[0], self.prefixes[1] + self.types[1], '1']
        if self.nbs:
            mod_list.extend([sigma + self.dpar, epsilon])
        if self.nbe:
            mod_list.extend([sigma, epsilon + self.dpar])
        returnlist = [fstring.format(*cont_list)] if add_original else []
        returnlist.append(fstring.format(*mod_list))
        return returnlist

    def prefix_type(self, top, tnames, tres, tresnames, prefix_number):
        """
        looks for a line in topology section 'atoms'
        and modifies it as needed
        :return: None
        """
        for x in range(len(tnames[prefix_number])):
            name = tnames[prefix_number][x]
            resname = tresnames[prefix_number][x]
            resnum = tres[prefix_number][x]
            for section in top.list_sections('atoms'):
                for line in range(len(top.sections[section])):
                    lspl = top.sections[section][line].split()
                    if len(lspl) >= 7 and lspl[4] == name and int(lspl[2]) == int(resnum) \
                            and lspl[3] == resname and lspl[1] == self.types[prefix_number]:
                        modline = self.add_prefix_to_line(lspl, prefix_number)
                        top.sections[section][line] = modline

    def add_prefix_to_line(self, cont_list, prefix_number):
        """
        takes a line from the topology section 'atoms'
        and increases charge by self.dpar;
        also changes atomtype to modified (Y-prefixed)
        :param cont_list: line passed as a list of non-whitespace strings
        :param prefix_number: int, which prefix to choose
        :return: an assembled line with type and/or charge modified
        """
        prefix = self.prefixes[prefix_number]
        cont_list[6] = float(cont_list[6])
        if not self.alch:
            cont_list[1] = prefix + cont_list[1]
            fstring = '{:>6s}{:>11s}{:>7s}{:>7s}{:>7s}{:>7s}{:>11.4f}{:>11s}\n'
            if len(cont_list) == 7:
                cont_list = cont_list + [' ']
            return fstring.format(*cont_list[:8])
        else:
            fstring = '{:>6s}{:>11s}{:>7s}{:>7s}{:>7s}{:>7s}{:>11.4f}{:>11s}{:>11s}{:>11.4f}{:>11s}\n'
            if len(cont_list) == 7:
                raise RuntimeError("You need to provide mass for the atom to be changed")
            return fstring.format(*cont_list[:8], prefix + cont_list[1], cont_list[6], cont_list[7])

    def types_are_consistent(self, top, tnames, tres, tresnames):
        """
        if many atoms are passed in a single input line, need to make sure
        they are all of the same type; raises ValueError if this is not the case
        :return: None
        """
        consistent1 = all([self.get_type(top, tnames[0][0], tres[0][0], tresnames[0][0])
                           == self.get_type(top, tnames[0][n], tres[0][n], tresnames[0][n])
                           for n in range(len(tnames[0]))])
        consistent2 = all([self.get_type(top, tnames[1][0], tres[1][0], tresnames[1][0])
                           == self.get_type(top, tnames[1][n], tres[1][n], tresnames[1][n])
                           for n in range(len(tnames[1]))])
        if not consistent1 or not consistent2:
            raise ValueError("atoms within a single nbfix modification need to have consistent types")

    @staticmethod
    def add_nbfix_section(top):
        position = top.section_headers.index('moleculetype')
        top.section_headers.insert(position, 'nonbond_params')
        top.sections.insert(position, ['; NBFIX  rmin=<charmm_rmin>/2^(1/6), eps=4.184*<charmm_eps>\n',
                                       ';name   type1  type2  1  sigma   epsilon\n', '[ nonbond_params ]\n'])

    def check_lj_params(self, top):
        for section in top.list_sections('atomtypes'):
            for line in range(len(top.sections[section])):
                lspl = top.sections[section][line].split()
                if len(lspl) > 6 and lspl[0] in self.types:
                    self.LJ_sigmas.append(float(lspl[5]))
                    self.LJ_epsilons.append(float(lspl[6]))
        if len(self.LJ_sigmas) == 2 and len(self.LJ_epsilons) == 2:
            pass
        elif len(self.LJ_sigmas) == 1 and len(self.LJ_epsilons) == 1 and self.types[0] == self.types[1]:
            self.LJ_sigmas.append(self.LJ_sigmas[0])
            self.LJ_epsilons.append(self.LJ_epsilons[0])
        else:
            raise RuntimeError("Found {} LJ sigmas and {} LJ epsilons for types {}, "
                               "check your topology".format(len(self.LJ_sigmas), len(self.LJ_sigmas), self.types))
        for section in top.list_sections('nonbond_params'):
            for line in range(len(top.sections[section])):
                lspl = top.sections[section][line].split()
                if len(lspl) > 4 and ((lspl[0], lspl[1]) == (self.types[0], self.types[1]) or
                                      (lspl[0], lspl[1]) == (self.types[1], self.types[0])):
                    self.NB_sigmas = float(lspl[3])
                    self.NB_epsilons = float(lspl[4])


class ModParam(Mod):
    def __init__(self, top, master, resnr='', resnames='', names='', changes=''):
        """
        Subclass to work with modified bonded params (angles and dihedrals).
        The class allows a lists of tuples of atoms that each corresponds
        to a single bonded term, angle or dihedral (all need to have
        the same sets of atomtypes for the calculation to be meaningful);
        the corresponding parameters will be modified in the molecule definition.
        The modification can be applied to all atoms bearing the specific type
        as well as only a subgroup of these.
        :param top: Top, a Top object corresponding to the modified topology
        :param resnr: list, a list of tuples, each containing a list of 3 or 4 residue numbers
        :param resnames: list, a list of tuples, each containing a list of 3 or 4 residue names
        :param names: list, a list of tuples, each containing a list of 3 or 4 atom names
        :param changes: str, single character denoting the change to be performed
        """
        super(ModParam, self).__init__(top, master, resnr, resnames, names)
        self.dih = 'd' if 'd' in changes else ''  # implement add to top line
        self.ang = 'a' if 'a' in changes else ''  # implement add to top line
        if self.dih and self.ang:
            raise RuntimeError("Please use separate Mods entries for angles and dihedrals")
        if self.dih and not self.ang:
            if len(changes) == 2:
                self.period = changes[1]
            else:
                self.period = ''
        self.npar = 3 if self.ang else 4
        if not (self.dih or self.ang):
            raise ValueError("no mode selected appropriate for an ModParam object")
        for t in self.top:
            self.types_are_consistent(t)
            self.types = tuple(self.get_type(t, n, r, rn) for n, r, rn in
                               zip(self.names[0], self.res[0], self.resnames[0]))
            self.origpar = []
            self.is_wildcard = False
            if self.ang:
                self.find_ang(t)
                self.mod_ang(t)
            if self.dih:
                self.find_dih(t)
                self.mod_dih(t)

    def mod_ang(self, top):
        """
        looks for lines in topology section 'angles'
        and modifies them as needed
        :return: None
        """
        sections = top.list_sections('angles')
        assert len(self.types) == 3
        for s in sections:
            molname = [x for x in top.mols.keys() if s in top.mols[x]][0]
            for l in range(len(top.sections[s])):
                lspl = top.sections[s][l].split()
                if len(lspl) > 3 and lspl[0][0] not in '[;':
                    try:
                        line_nums = [tuple(top.names_to_nums[molname]["{}-{}-{}".format(xr, xrn, xn)]
                                           for xr, xrn, xn in zip(r, rn, n))
                                     for r, rn, n in zip(self.resnames, self.res, self.names)]
                        if any([self.check_line(lspl[:3], line) for line in line_nums]):
                            modline = self.increment_angle(lspl)
                            top.sections[s][l] = modline
                    except KeyError:
                        pass

    def increment_angle(self, cont_list):
        """
        adds the modified angle parameters to specified line
        from topology section 'angles'
        :param cont_list: list, line to be modified, passed as a list of non-whitespace strings
        :return: str, modified and formatted line
        """
        assert len(self.origpar) == 1
        cont_list.extend(self.origpar[0]) # TODO fix
        fstring = '{:>5s} {:>5s} {:>5s} {:>5s}' + len(self.origpar[0]) * '{:>14} '
        return fstring.format(*cont_list)

    def mod_dih(self, top):
        """
        looks for lines in topology section 'dihedrals'
        and modifies them as needed
        :return: None
        """
        sections = top.list_sections('dihedrals')
        assert len(self.types) == 4
        for s in sections:
            molname = [x for x in top.mols.keys() if s in top.mols[x]][0]
            for l in range(len(top.sections[s])):
                lspl = top.sections[s][l].split()
                if len(lspl) > 4 and lspl[0][0] not in '[;':
                    try:
                        line_nums = [tuple(top.names_to_nums[molname]["{}-{}-{}".format(xr, xrn, xn)]
                                           for xr, xrn, xn in zip(r, rn, n))
                                     for r, rn, n in zip(self.resnames, self.res, self.names)]
                        if any([self.check_line(lspl[:4], line) for line in line_nums]):
                            modline = self.increment_dihedral(lspl)
                            top.sections[s][l] = modline
                    except KeyError:
                        pass

    def increment_dihedral(self, cont_list):
        """
        adds the modified dihedral parameters to specified line
        from topology section 'dihedrals'
        :param cont_list: list, line to be modified, passed as a list of non-whitespace strings
        :return: str, modified and formatted line
        """
        assert self.origpar
        result = ''
        fstring = '{:>5s} {:>5s} {:>5s} {:>5s} {:>5s}' + 6 * '{:>14} ' + '\n'
        for par in self.origpar:
            cont_temp = cont_list[:]
            opar = par[:]
            if (self.period and par[2] == self.period) or not self.period:
                par[1] = float(par[1]) + self.dpar
            cont_temp.extend(opar[:3] + par[:3])
            result += fstring.format(*cont_temp)
        return result

    def check_line(self, list_line_to_check, list_entries):
        """
        checks if two lines share the same parameter entries
        (possibly in reverse order, and including wildcards)
        :param list_line_to_check: list, line from topology file split at whitespaces
        :param list_entries: list, contains desired atomtypes as strings
        :return: bool, True if match found
        """
        n = len(list_entries)
        ok = all([self.eq_x(list_line_to_check[i], list_entries[i]) for i in range(n)]) \
            or all([self.eq_x(list_line_to_check[i], list_entries[n - i - 1]) for i in range(n)])
        return ok

    def find_ang(self, top):
        """
        looks for the desired angle parameters in the topology
        section "angletypes"
        :return: None
        """
        section = top.list_sections('angletypes')[0]
        for l in range(len(top.sections[section])):
            lspl = top.sections[section][l].split()
            if len(lspl) > 3 and lspl[0][0] not in '[;' and self.check_line(lspl, self.types):
                self.origpar.append(lspl[4:])
                self.origpar[-1][1] = float(self.origpar[-1][1]) + self.dpar

    def find_dih(self, top):
        """
        looks for the desired dihedral parameters in the topology
        section "dihedraltypes"; also accounts for cases in which
        wildcarded and non-wildcard params can match the request,
        overwriting wildcards with non-wildcards whenever possible
        and saving wildcards when no alternatives exist
        :return: None
        """
        sections = top.list_sections('dihedraltypes')
        for s in sections:
            for l in range(len(top.sections[s])):
                lspl = top.sections[s][l].split()
                if len(lspl) > 6 and lspl[0][0] not in '[;' and self.check_line(lspl, self.types) \
                        and (not self.period or lspl[7] == self.period):
                    if ('X' in lspl[:4] and not self.origpar) or self.is_wildcard:
                        self.origpar.append(lspl[5:])
                        self.is_wildcard = True
                    else:
                        self.origpar.append(lspl[5:])
        if not self.period:
            if not len(self.origpar) == 1:
                raise RuntimeError('Found multiple dihedral entries for the requested modification. To work with '
                                   'multiple dihedral entries, use notation "d1", "d2" etc. where 1, 2, ... denotes'
                                   ' the periodicity of the given term.')

    def types_are_consistent(self, top):
        """
        checks whether all input atoms share the same type signature,
        (t1 t2 t3 t4) for dihedrals or (t1 t2 t3) for angles;
        raises ValueError if not the case
        :return: None
        """
        types_list = [tuple(self.get_type(top, qn, qrn, qr) for qn, qrn, qr in zip(n, rn, r))
                      for n, rn, r in zip(self.names, self.res, self.resnames)]
        ok = all([self.check_line(types_list[0], types_list[n]) for n in range(1, len(types_list))])
        if not ok:
            raise ValueError("atoms within a single parameter modification need to have consistent types")

    @staticmethod
    def eq_x(one, another):
        """
        compares the identity of two types, allowing for wildcards
        ("X" can match any type in the "dihedraltypes" section)
        :param one: str, first element to compare
        :param another: str, second element to compare
        :return: bool, True if match
        """
        if one == another or one == "X" or another == "X":
            return True
        else:
            return False
